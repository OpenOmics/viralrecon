#!/usr/bin/env python3
"""Stage driver for the POD5-only nf-core/viralrecon TCR extension.

External tools are invoked with argument lists, never interpolated shell commands.
BAM is authoritative throughout, retaining Dorado read-group/model/move metadata.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
from collections import Counter, defaultdict

from tcr_pod5_targets import Endpoint, infer_targets, load_intervals, target_candidate

from tcr_pod5_core import (VERSION, Settings, Sam, assign, build_scheme, evidence,
                            full_primer_match, rc, read_fasta, screen, write_fasta)


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def run(cmd, output=None, log=None):
    cmd = [str(x) for x in cmd]
    print("COMMAND " + json.dumps(cmd), file=sys.stderr)
    with (open(output, "wb") if output else open(os.devnull, "wb")) as out:
        if log:
            with open(log, "ab") as err:
                subprocess.run(cmd, stdout=out, stderr=err, check=True)
        else:
            subprocess.run(cmd, stdout=out, check=True)


def capture(cmd):
    return subprocess.check_output([str(x) for x in cmd], text=True)


def dorado_version(executable, expected):
    p = subprocess.run([executable, "--version"], text=True, capture_output=True, check=True)
    text = p.stdout + p.stderr
    m = re.search(r"(?<!\d)(\d+\.\d+\.\d+)", text)
    if not m or m.group(1) != expected:
        raise ValueError(f"Require Dorado {expected}; received {text.strip()!r}")
    return text.strip()


def sam_lines(path, samtools):
    if str(path).endswith(".sam"):
        with open(path) as f:
            yield from f
        return
    p = subprocess.Popen([samtools, "view", "--no-PG", "-h", str(path)], stdout=subprocess.PIPE, text=True)
    try:
        yield from p.stdout
    finally:
        p.stdout.close()
        result = p.wait()
        if result:
            raise RuntimeError(f"samtools view failed for {path} (exit {result})")


def alignments(path, samtools):
    for line in sam_lines(path, samtools):
        if not line.startswith("@"):
            yield Sam.parse(line)


def sam_to_bam(path, output, samtools):
    run([samtools, "view", "-b", "-o", output, path])


def bam_fastq(path, output, samtools):
    # Input is unaligned single-end BAM. Preserve the read UUID as the first token.
    with gzip.open(output, "wt") as f:
        for r in alignments(path, samtools):
            if r.flag & (256 | 2048):
                continue
            if r.seq == "*" or r.qual == "*":
                raise ValueError(f"Missing basecalls or qualities for {r.name}")
            if r.flag & 16:
                seq, qual = rc(r.seq), r.qual[::-1]
            else:
                seq, qual = r.seq, r.qual
            tag = " ".join(f"{k}={r.tags[k]}" for k in ("BC", "RG", "TS") if k in r.tags)
            f.write(f"@{r.name} {tag}\n{seq}\n+\n{qual}\n")


def cmd_audit(a):
    import pod5  # Deliberately only required for the POD5 audit stage.
    entries = json.loads(Path(a.manifest).read_text())
    db = sqlite3.connect("read_ids.sqlite")
    db.execute("CREATE TABLE ids (uuid TEXT PRIMARY KEY, batch TEXT, source TEXT)")
    rows, duplicates, empty, batch_counts = [], [], [], Counter()
    with open("pod5_read_ids.tsv", "w") as ids:
        ids.write("read_id\tbatch\tacquisition_id\tflow_cell_id\tsample_id\tsource\n")
        for entry in entries:
            path, batch = Path(entry["staged"]), entry["batch"]
            source = entry["source"]
            if path.stat().st_size == 0:
                empty.append(dict(batch=batch, source=source))
                continue
            metadata, count = set(), 0
            # POD5 container and read-table validation; signal decoding happens in Dorado.
            with pod5.Reader(path) as reader:
                for r in reader.reads():
                    rid = str(r.read_id)
                    info = r.run_info
                    key = tuple(str(getattr(info, k, "")) for k in
                                ("acquisition_id", "flow_cell_id", "sample_id", "sequencing_kit", "sample_rate"))
                    metadata.add(key)
                    if entry.get("kit") and key[3] and entry["kit"].upper() != key[3].upper():
                        raise ValueError(f"Sequencing-kit mismatch for {source}: manifest={entry['kit']}, POD5={key[3]}")
                    try:
                        db.execute("INSERT INTO ids VALUES (?,?,?)", (rid, batch, source))
                    except sqlite3.IntegrityError:
                        old = db.execute("SELECT batch,source FROM ids WHERE uuid=?", (rid,)).fetchone()
                        duplicates.append(dict(read_id=rid, first_batch=old[0], first_source=old[1],
                                               second_batch=batch, second_source=source))
                    ids.write("\t".join((rid, batch, *key[:3], source)) + "\n")
                    count += 1
            if count == 0:
                empty.append(dict(batch=batch, source=source))
            batch_counts[batch] += count
            rows.append(dict(batch=batch, source=source, bytes=path.stat().st_size,
                             reads=count, run_metadata=sorted(metadata)))
        db.commit()
    db.close()
    dump("pod5_audit.json", dict(files=rows, empty_files=empty, duplicate_reads=duplicates,
                                  batch_reads=dict(batch_counts), pod5_version=pod5.__version__))
    if duplicates:
        raise ValueError(f"Found {len(duplicates)} repeated POD5 read UUIDs. Resolve provenance; do not silently deduplicate.")
    if empty and not a.allow_empty:
        raise ValueError(f"Found {len(empty)} empty POD5 files. See pod5_audit.json. Explicit --tcr_allow_empty_pod5 is required to ignore them.")
    for batch in {x["batch"] for x in entries}:
        if not batch_counts[batch]:
            raise ValueError(f"Batch {batch} contains no POD5 reads")
        rates = {m[4] for row in rows if row["batch"] == batch for m in row["run_metadata"]}
        if len(rates) > 1:
            raise ValueError(f"Batch {batch} has multiple signal sample rates; separate by chemistry/model: {rates}")
    Path("audit.ok").write_text("POD5 read metadata validated; no duplicate read UUIDs.\n")


def cmd_prepare(a):
    expected = [x.strip() for x in a.expected_constructs.split(";")] if a.expected_constructs else None
    pairs = build_scheme(a.primers, a.references, a.outdir,
                         a.expected_pairs, a.barcode_errors, expected)
    out = Path(a.outdir)
    for name, file in (("primers", a.primers), ("references", a.references)):
        shutil.copy2(file, out / ("input_" + name + Path(file).suffix))
    dump(out / "input_checksums.json", {
        str(file): hashlib.sha256(Path(file).read_bytes()).hexdigest()
        for file in (a.primers, a.references)
    })
    print(f"Prepared {len(pairs)} legal well pairs", file=sys.stderr)


def cmd_basecall(a):
    version = dorado_version(a.dorado, a.version)
    files = json.loads(Path(a.files).read_text())
    model = Path(a.model)
    if not model.is_dir() or not (model / "config.toml").is_file():
        raise ValueError("Provide a pre-downloaded, versioned Dorado simplex model directory containing config.toml")
    staging = Path("valid_pod5")
    staging.mkdir()
    for i, p in enumerate(files):
        p = Path(p)
        if p.stat().st_size:
            (staging / f"part_{i:06d}.pod5").symlink_to(p.resolve())
    if not list(staging.iterdir()):
        raise ValueError("No nonempty POD5 files to basecall")
    cmd = [a.dorado, "basecaller", model, staging, "--device", a.device,
           "--no-trim", "--emit-moves", "--min-qscore", "0", "--disable-read-splitting"]
    run(cmd, "calls.bam", "dorado_basecall.log")
    run([a.samtools, "quickcheck", "-u", "calls.bam"])
    count = int(capture([a.samtools, "view", "-c", "calls.bam"]).strip())
    expected = json.loads(Path(a.audit).read_text())["batch_reads"][a.batch]
    if count != expected:
        raise ValueError(f"POD5/basecall read-count mismatch: {expected} vs {count}; see Dorado log")
    run([a.dorado, "summary", "calls.bam"], "basecall_summary.tsv", "dorado_summary.log")
    dump("basecall_provenance.json", dict(batch=a.batch, command=[str(x) for x in cmd],
          dorado_version=version, reads=count, model=str(model.resolve()),
          model_config_sha256=hashlib.sha256((model / "config.toml").read_bytes()).hexdigest(),
          read_splitting=False, basecall_quality_filter=0))


def barcode_well(label, n_to_well):
    if not label or label == "unclassified":
        return "unclassified"
    # Known Dorado custom-kit label encodings. Never infer identity from chunk filenames.
    m = re.fullmatch(r"(?:(?:TCR_WELLS|TCR)_)?(?:barcode|FBC|RBC)(\d+)", label)
    if not m or int(m.group(1)) not in n_to_well:
        raise ValueError(f"Unrecognized custom barcode BC tag {label!r}; inspect demux output before changing mapping")
    return n_to_well[int(m.group(1))]


def cmd_demux(a):
    dorado_version(a.dorado, a.version)
    scheme = Path(a.scheme)
    pairs = json.loads((scheme / "scheme.json").read_text())["pairs"]
    n_to_well = {p["barcode_number"]: p["well"] for p in pairs}
    run([a.dorado, "demux", a.bam, "--kit-name", "TCR_WELLS",
         "--barcode-arrangement", scheme / "barcodes.toml",
         "--barcode-sequences", scheme / "barcodes.fasta",
         "--barcode-both-ends", "--no-trim", "--emit-summary",
         "--threads", a.threads, "--output-dir", "raw_demux"], log="dorado_demux.log")
    bams = sorted(Path("raw_demux").rglob("*.bam"))
    if not bams:
        raise ValueError("Dorado demux produced no BAM files")
    headers = {}
    for b in bams:
        for line in capture([a.samtools, "view", "--no-PG", "-H", b]).splitlines(keepends=True):
            if line.startswith("@HD"):
                continue
            idmatch = re.search(r"(?:^|\t)ID:([^\t\n]+)", line)
            key = (line[:3], idmatch.group(1) if idmatch else line)
            if key in headers and headers[key] != line:
                raise ValueError(f"Conflicting SAM header for {key}")
            headers[key] = line
    Path("demux").mkdir()
    Path("demux/wells").mkdir()
    names = [p["well"] for p in pairs] + ["unclassified"]
    handles = {w: open(f"demux/{w}.sam", "w") for w in names}
    for f in handles.values():
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        f.writelines(headers.values())
    counts, seen = Counter(), set()
    with open("demux/read_assignments.tsv", "w") as assignments:
        assignments.write("read_id\tbatch\twell\tBC\n")
        try:
            for b in bams:
                for r in alignments(b, a.samtools):
                    if r.name in seen:
                        raise ValueError(f"Duplicate demultiplexed read UUID: {r.name}")
                    seen.add(r.name)
                    label = r.tags.get("BC", "unclassified")
                    well = barcode_well(label, n_to_well)
                    counts[well] += 1
                    handles[well].write(r.line)
                    assignments.write(f"{r.name}\t{a.batch}\t{well}\t{label}\n")
        finally:
            for f in handles.values():
                f.close()
    expected = int(capture([a.samtools, "view", "-c", a.bam]).strip())
    if len(seen) != expected:
        raise ValueError(f"Read accounting failed during demultiplexing: {expected} input, {len(seen)} output")
    rows = []
    for well in names:
        sam = Path(f"demux/{well}.sam")
        bam = Path(f"demux/wells/{well}.bam")
        sam_to_bam(sam, bam, a.samtools)
        sam.unlink()
        rows.append(dict(well=well, reads=counts[well], batch=a.batch))
    bam_fastq("demux/wells/unclassified.bam", "demux/unclassified.fastq.gz", a.samtools)
    dump("demux/manifest.json", rows)
    shutil.copy2("dorado_demux.log", "demux/dorado_demux.log")
    for p in Path("raw_demux").rglob("*.txt"):
        shutil.copy2(p, Path("demux") / ("dorado_" + p.name))
    for p in Path("raw_demux").rglob("*.tsv"):
        shutil.copy2(p, Path("demux") / ("dorado_" + p.name))
    dump("demux/settings.json", dict(dorado_version=a.version, both_ends=True, trim=False))


def check_trim(original: Sam, trimmed: Sam, fwd: str, rev: str, cfg: Settings):
    if original.seq == "*" or trimmed.seq == "*" or not trimmed.seq:
        return "missing_sequence", None
    hits = []
    for direction, left, right in (("+", fwd, rc(rev)), ("-", rev, rc(fwd))):
        le = full_primer_match(left, original.seq[:cfg.primer_window])
        offset = max(0, len(original.seq) - cfg.primer_window)
        ri = full_primer_match(right, original.seq[offset:])
        if le[0] <= cfg.primer_errors and ri[0] <= cfg.primer_errors:
            hits.append((le[0] + ri[0], direction, le[2], ri[1] + offset))
    if not hits:
        return "two_full_primers_not_confirmed", None
    hits.sort()
    if len(hits) > 1 and hits[0][0] == hits[1][0]:
        return "ambiguous_orientation", None
    _, direction, left_end, right_start = hits[0]
    start = original.seq.find(trimmed.seq)
    if start < 0 or original.seq.find(trimmed.seq, start + 1) >= 0:
        return "trimmed_sequence_not_unique_substring", direction
    end = start + len(trimmed.seq)
    if not (left_end <= start <= left_end + cfg.max_overtrim and
            right_start - cfg.max_overtrim <= end <= right_start):
        return "primer_boundaries_not_removed_as_expected", direction
    if original.qual[start:end] != trimmed.qual:
        return "trimmed_quality_mismatch", direction
    return "pass", direction


def cmd_trim(a):
    cfg = Settings.from_json(a.settings)
    dorado_version(a.dorado, a.version)
    primers = read_fasta(a.primers)
    fwd, rev = primers[a.well + "_front"], primers[a.well + "_rear"]
    run([a.dorado, "trim", a.bam, "--sequencing-kit", a.kit,
         "--primer-sequences", a.primers, "--threads", a.threads],
        "dorado_trimmed.bam", "dorado_trim.log")
    raw = {r.name: r for r in alignments(a.bam, a.samtools)}
    seen, counts = set(), Counter()
    with open("clean.sam", "w") as good, open("trim_failed.sam", "w") as bad, open("trim_audit.tsv", "w") as tab:
        tab.write("read_id\tstatus\torientation\tbefore_length\tafter_length\n")
        for line in sam_lines("dorado_trimmed.bam", a.samtools):
            if line.startswith("@"):
                good.write(line)
                bad.write(line)
                continue
            r = Sam.parse(line)
            if r.name in seen or r.name not in raw:
                raise ValueError(f"Unexpected or repeated UUID after trimming: {r.name}")
            seen.add(r.name)
            result, direction = check_trim(raw[r.name], r, fwd, rev, cfg)
            counts[result] += 1
            (good if result == "pass" else bad).write(line)
            tab.write(f"{r.name}\t{result}\t{direction or 'unknown'}\t{len(raw[r.name].seq)}\t{len(r.seq)}\n")
    if seen != set(raw):
        raise ValueError("Dorado trim lost records; inspect command output")
    sam_to_bam("clean.sam", "clean.bam", a.samtools)
    sam_to_bam("trim_failed.sam", "trim_failed.bam", a.samtools)
    bam_fastq("clean.bam", f"{a.well}.fastq.gz", a.samtools)
    dump("trim_stats.json", dict(input_reads=len(raw), passed_reads=counts["pass"],
                                 failed_reads=len(raw) - counts["pass"], reasons=dict(counts)))
    Path("clean.sam").unlink()
    Path("trim_failed.sam").unlink()


def do_align(dorado, samtools, refs, bam, prefix, threads, secondary=True):
    options = "-x lr:hq --eqx --secondary " + ("yes -N 50" if secondary else "no")
    count = int(capture([samtools, "view", "-c", bam]).strip())
    if count:
        run([dorado, "aligner", refs, bam, "--threads", threads, "--mm2-opts", options],
            prefix + ".unsorted.bam", prefix + ".align.log")
    else:
        # Empty wells still yield valid alignment headers and an explicit QC decision.
        with open(prefix + ".sam", "w") as f:
            f.write("@HD\tVN:1.6\tSO:unsorted\n")
            for name, seq in read_fasta(refs).items():
                f.write(f"@SQ\tSN:{name}\tLN:{len(seq)}\n")
        sam_to_bam(prefix + ".sam", prefix + ".unsorted.bam", samtools)
        Path(prefix + ".sam").unlink()
    run([samtools, "sort", "-@", threads, "-o", prefix + ".bam", prefix + ".unsorted.bam"])
    run([samtools, "index", prefix + ".bam"])
    Path(prefix + ".unsorted.bam").unlink()


def cmd_align(a):
    dorado_version(a.dorado, a.version)
    do_align(a.dorado, a.samtools, a.references, a.bam, "aligned", a.threads)
    run([a.samtools, "sort", "-n", "-@", a.threads, "-o", "name_sorted.bam", "aligned.bam"])


def cmd_infer_targets(a):
    cfg = Settings.from_json(a.settings)
    refs = read_fasta(a.references)
    manifest = json.loads(Path(a.manifest).read_text())
    if manifest["batch"] != a.batch:
        raise ValueError("Inference manifest batch mismatch")
    seen, wells, points, excluded, assigned = set(), set(), [], Counter(), Counter()
    for item in sorted(manifest["alignments"], key=lambda item: item["well"]):
        well = item["well"]
        if well in wells:
            raise ValueError(f"Repeated well in target inference: {well}")
        wells.add(well)
        for name, records in itertools.groupby(alignments(item["bam"], a.samtools), key=lambda r: r.name):
            if name in seen:
                raise ValueError(f"Repeated read UUID across batch alignments: {name}")
            seen.add(name)
            status, chosen = assign(list(records), refs, cfg, require_target_coverage=False)
            if chosen is None:
                excluded[status] += 1
                continue
            assigned[status] += 1
            reason = target_candidate(chosen, cfg)
            if reason:
                excluded[reason] += 1
                continue
            points.append(Endpoint(name, well, chosen.ref, chosen.start, chosen.end, chosen.strand))
    report, used = infer_targets(refs, points, cfg, a.batch, manifest.get("expected_constructs"))
    report.update(input_reads=len(seen), input_wells=len(wells),
                  competitive_assignments=dict(assigned), excluded_reads=dict(excluded))
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=False)
    dump(out / "inferred_targets.json", report)
    with open(out / "inferred_targets.bed", "w") as bed, open(out / "target_summary.tsv", "w") as table:
        fields = ["reference", "status", "consensus_eligible", "start", "end", "candidate_reads",
                  "candidate_wells", "supporting_reads", "supporting_wells", "balanced_cluster_fraction", "reason"]
        writer = csv.DictWriter(table, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for name, row in sorted(report["references"].items()):
            writer.writerow(row)
            if row["status"] == "INFERRED":
                bed.write(f"{name}\t{row['start']}\t{row['end']}\t{name}\n")
    with open(out / "boundary_evidence.tsv", "w") as table:
        table.write("read_id\twell\treference\tstart\tend\tstrand\tsupports_inferred_target\n")
        for row in sorted(points, key=lambda row: (row.reference, row.well, row.read_id)):
            table.write(f"{row.read_id}\t{row.well}\t{row.reference}\t{row.start}\t{row.end}\t{row.strand}\t{int(row.read_id in used)}\n")


def cmd_qc(a):
    cfg = Settings.from_json(a.settings)
    refs = read_fasta(a.references)
    target_report = json.loads(Path(a.inferred_targets).read_text())
    targets = load_intervals(target_report, refs, a.batch)
    counts, good, seen = Counter(), [], set()
    Path("qc").mkdir()
    with open("qc/read_assignments.tsv", "w") as table:
        table.write("read_id\tassignment\tstrand\n")
        for name, records in itertools.groupby(alignments(a.bam, a.samtools), key=lambda r: r.name):
            if name in seen:
                raise ValueError("QC requires a query-name sorted BAM; repeated nonconsecutive read ID")
            seen.add(name)
            status, chosen = assign(list(records), refs, cfg, require_target_coverage=False)
            counts[status] += 1
            if chosen:
                good.append(chosen)
            table.write(f"{name}\t{status}\t{chosen.strand if chosen else '.'}\n")
    trim = json.loads(Path(a.trim_stats).read_text())
    if sum(counts.values()) != trim["passed_reads"]:
        raise ValueError("Alignment/trimmed read accounting mismatch")
    report, ev, selected = screen(list(counts.elements()), good, refs, cfg, trim,
                                  targets=targets, expected_constructs=target_report["expected_constructs"])
    report.update(sample=a.sample, batch=a.batch, well=a.well, trim_stats=trim,
                  target_inference=target_report["references"][report["dominant_construct"]])
    dump("qc/decision.json", report)
    with open("qc/allele_support.tsv", "w") as table:
        table.write("kind\tposition_0based\thomopolymer\tallele_counts\n")
        for site in sorted(ev["sites"], key=lambda s: (s["pos"], s["kind"])):
            table.write(f"{site['kind']}\t{site['pos']}\t{site['homopolymer']}\t{json.dumps(site['alleles'], sort_keys=True)}\n")
    with open("qc/candidate_haplotypes.tsv", "w") as table:
        table.write("read_id\tstrand\tobserved_candidate_alleles\n")
        for r in selected:
            data = {f"{s['kind']}:{s['pos']}": s["calls"].get(r.name, ".") for s in ev["sites"]}
            table.write(f"{r.name}\t{r.strand}\t{json.dumps(data, sort_keys=True)}\n")
    if report["status"] == "PASS_SCREEN":
        Path("pass").mkdir()
        write_fasta("pass/draft.fasta", {a.sample: ev["draft"]})
        Path("pass/accepted.ids").write_text("".join(r.name + "\n" for r in selected))
        shutil.copy2("qc/decision.json", "pass/decision.json")


def check_basecall_models(path, samtools):
    models = set()
    for line in capture([samtools, "view", "--no-PG", "-H", path]).splitlines(keepends=True):
        if line.startswith("@RG"):
            m = re.search(r"basecall_model=([^\s]+)", line)
            if not m:
                raise ValueError("Read-group header lacks a basecall_model; refusing --ignore-read-groups")
            models.add(m.group(1))
    if len(models) != 1:
        raise ValueError(f"Polishing requires exactly one basecall model across all RGs; found {models}")
    return next(iter(models))


def cmd_polish(a):
    cfg = Settings.from_json(a.settings)
    dorado_version(a.dorado, a.version)
    passed = Path(a.passed)
    if json.loads((passed / "decision.json").read_text())["status"] != "PASS_SCREEN":
        raise ValueError("Refusing to make consensus for a non-passing well")
    run([a.samtools, "view", "-b", "-N", passed / "accepted.ids", "-o", "accepted.bam", a.bam])
    model = check_basecall_models("accepted.bam", a.samtools)
    do_align(a.dorado, a.samtools, passed / "draft.fasta", "accepted.bam", "draft_aligned", a.threads, False)
    run([a.dorado, "polish", "draft_aligned.bam", passed / "draft.fasta",
         "--ignore-read-groups", "--device", a.device, "--threads", a.threads,
         "--infer-threads", a.threads, "--models-directory", a.models],
        "polished.fasta", "dorado_polish.log")
    do_align(a.dorado, a.samtools, "polished.fasta", "accepted.bam", "polished_aligned", a.threads, False)
    refs = read_fasta("polished.fasta")
    if len(refs) != 1:
        raise ValueError("Expected exactly one polished amplicon per passing well")
    name, seq = next(iter(refs.items()))
    groups = defaultdict(list)
    for r in alignments("polished_aligned.bam", a.samtools):
        groups[r.name].append(r)
    reads, post_status = [], Counter()
    for records in groups.values():
        status, chosen = assign(records, refs, cfg)
        post_status[status] += 1
        if chosen:
            reads.append(chosen)
    ev = evidence(reads, seq, cfg)
    output, masked = [], []
    for p, base in enumerate(seq):
        counts = Counter(ev["pile"].get(p, {}).values())
        depth = sum(counts.values())
        if base in "ACGT" and depth >= cfg.min_depth and counts[base] / max(1, depth) >= cfg.consensus_fraction:
            output.append(base)
        else:
            output.append("N")
            masked.append(p)
    fixed_indels = []
    for event, calls in ev["events"]:
        counts = Counter(calls.values())
        if sum(counts.values()) >= cfg.min_depth and counts["ALT"] / sum(counts.values()) >= cfg.consensus_fraction:
            fixed_indels.append(event)
    reasons = []
    unresolved = sum(post_status.values()) - len(reads)
    if unresolved / max(1, sum(post_status.values())) > cfg.max_unresolved_fraction:
        reasons.append("REVIEW_POSTPOLISH_ALIGNMENTS")
    if post_status.get("split_alignment", 0) >= cfg.alert_reads:
        reasons.append("REVIEW_POSTPOLISH_SPLIT_ALIGNMENTS")
    if ev["sites"]:
        reasons.append("REVIEW_POSTPOLISH_VARIATION")
    if fixed_indels:
        reasons.append("REVIEW_POSTPOLISH_INDEL_DISAGREEMENT")
    if 1 - len(masked) / len(seq) < cfg.min_callable_fraction:
        reasons.append("REVIEW_POSTPOLISH_COVERAGE")
    if len(reads) < cfg.min_reads:
        reasons.append("REVIEW_POSTPOLISH_READ_SUPPORT")
    Path("consensus_result").mkdir()
    with open("consensus_result/masked_positions.bed", "w") as f:
        for p in masked:
            f.write(f"{a.sample}\t{p}\t{p + 1}\n")
    report = dict(sample=a.sample, status="CONSENSUS_PASS" if not reasons else "REVIEW",
                  reasons=reasons, masked_bases=len(masked), length=len(seq), basecall_model=model,
                  residual_fixed_indels=fixed_indels, alignment_counts=dict(post_status), dorado_version=a.version,
                  coordinates="Masked BED uses the polished amplicon, not the original vector reference")
    dump("consensus_result/decision.json", report)
    shutil.copy2("dorado_polish.log", "consensus_result/dorado_polish.log")
    if not reasons:
        write_fasta("consensus_result/consensus.fasta", {a.sample: "".join(output)})


def cmd_report(a):
    data = json.loads(Path(a.manifest).read_text())
    rows, seqs = [], {}
    target_reports = [json.loads(Path(d, "inferred_targets.json").read_text()) for d in data.get("target_dirs", [])]
    if len({r["batch"] for r in target_reports}) != len(target_reports):
        raise ValueError("Duplicate batch target reports")
    qc = {json.loads(Path(p, "decision.json").read_text())["sample"]: p for p in data["qc_dirs"]}
    pol = {json.loads(Path(p, "decision.json").read_text())["sample"]: p for p in data["polish_dirs"]}
    for d in data["demux_dirs"]:
        for item in json.loads(Path(d, "manifest.json").read_text()):
            if item["well"] == "unclassified":
                continue
            sample = f"{item['batch']}__{item['well']}"
            row = dict(sample=sample, batch=item["batch"], well=item["well"], raw_assigned_reads=item["reads"],
                       status="NO_READS", reasons="No reads classified to this well", dominant_construct="",
                       usable_reads=0, callable_fraction="", target_status="", target_start="", target_end="")
            if item["reads"] and sample not in qc:
                raise ValueError(f"Missing well QC result: {sample}")
            if sample in qc:
                q = json.loads(Path(qc[sample], "decision.json").read_text())
                row.update(status=q["status"], reasons=";".join(q["reasons"]),
                           dominant_construct=q["dominant_construct"], usable_reads=q["assigned_reads"],
                           callable_fraction=q["callable_fraction"],
                           target_status=q.get("target_inference", {}).get("status", ""),
                           target_start=q.get("target_inference", {}).get("start"),
                           target_end=q.get("target_inference", {}).get("end"))
                if q["status"] == "PASS_SCREEN" and sample not in pol:
                    raise ValueError(f"Missing consensus result for passing well: {sample}")
            if sample in pol:
                p = json.loads(Path(pol[sample], "decision.json").read_text())
                row.update(status=p["status"], reasons=";".join(p["reasons"]))
                if p["status"] == "CONSENSUS_PASS":
                    s = read_fasta(Path(pol[sample], "consensus.fasta"))
                    if set(s) != {sample} or sample in seqs:
                        raise ValueError("Consensus identifier collision")
                    seqs.update(s)
            rows.append(row)
    rows.sort(key=lambda r: (r["batch"], r["well"][0], int(r["well"][1:])))
    with open("well_summary.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else [], delimiter="\t")
        if rows:
            w.writeheader()
            w.writerows(rows)
    write_fasta("accepted_consensus.fasta", dict(sorted(seqs.items())))
    dump("run_summary.json", dict(wells=len(rows), accepted_consensuses=len(seqs),
          statuses=dict(Counter(r["status"] for r in rows)), implementation_version=VERSION,
          batch_targets={r["batch"]: r for r in sorted(target_reports, key=lambda x: x["batch"])},
          note="Exploratory mixture screening. PASS does not establish one physical template."))
    print(f"Completed: {len(rows)} wells, {len(seqs)} accepted consensus sequences", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit")
    audit.add_argument("--manifest", required=True)
    audit.add_argument("--allow-empty", action="store_true")
    prep = sub.add_parser("prepare")
    for flag in ("primers", "references", "outdir"):
        prep.add_argument("--" + flag, required=True)
    prep.add_argument("--expected-constructs", default="",
                      help="Semicolon-separated FASTA IDs eligible for consensus; omitted means all")
    prep.add_argument("--expected-pairs", type=int, default=48)
    prep.add_argument("--barcode-errors", type=int, default=1)
    base = sub.add_parser("basecall")
    for flag in ("files", "model", "audit", "batch"):
        base.add_argument("--" + flag, required=True)
    base.add_argument("--device", default="cuda:0")
    demux = sub.add_parser("demux")
    for flag in ("bam", "scheme", "batch"):
        demux.add_argument("--" + flag, required=True)
    trim = sub.add_parser("trim")
    for flag in ("bam", "primers", "well", "kit", "settings"):
        trim.add_argument("--" + flag, required=True)
    aln = sub.add_parser("align")
    for flag in ("bam", "references"):
        aln.add_argument("--" + flag, required=True)
    inf = sub.add_parser("infer-targets")
    for flag in ("manifest", "references", "settings", "batch", "outdir"):
        inf.add_argument("--" + flag, required=True)
    qc = sub.add_parser("qc")
    qc.add_argument("--inferred-targets", required=True)
    for flag in ("bam", "references", "settings", "trim-stats", "sample", "batch", "well"):
        qc.add_argument("--" + flag, required=True)
    pol = sub.add_parser("polish")
    for flag in ("bam", "passed", "settings", "models", "sample"):
        pol.add_argument("--" + flag, required=True)
    pol.add_argument("--device", default="cpu")
    rep = sub.add_parser("report")
    rep.add_argument("--manifest", required=True)
    for tool in (base, demux, trim, aln, pol):
        tool.add_argument("--dorado", default="dorado")
        tool.add_argument("--version", default="2.1.2")
    for tool in (base, demux, trim, aln, qc, pol, inf):
        tool.add_argument("--samtools", default="samtools")
        tool.add_argument("--threads", default="4")
    a = p.parse_args()
    try:
        globals()["cmd_" + a.command.replace("-", "_")](a)
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as e:
        p.exit(1, f"ERROR: {e}\n")


if __name__ == "__main__":
    main()
