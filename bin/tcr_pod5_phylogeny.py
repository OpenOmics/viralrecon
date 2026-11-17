#!/usr/bin/env python3
"""Optional, batch-level MAFFT/RAxML-NG analysis of final TCR outputs.

No reads are recalled, no QC decisions are changed, and no reference bases are
inserted into a consensus. All coordinates in crop tables are 0-based half-open.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys

VERSION = "1.0.0"
DNA = set("ACGTRYSWKMBDHVN")
SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
LIMITATIONS = (
    "Exploratory sequence comparison, not a single-template assay or proof of natural ancestry. "
    "REVIEW sequences remain unaccepted, including masked unpolished drafts. "
    "N/ambiguous bases and alignment gaps do not establish reference identity; gaps are treated "
    "as missing data by the DNA model, not as a separate indel character. Engineered construct "
    "differences and shared vector sequence can dominate the tree. Identical tips are retained "
    "but their internal ordering is not evidence of a biological relationship. The tree is "
    "unrooted; references are comparison tips, not an assumed outgroup."
)


class NotAnalyzable(Exception):
    """Valid inputs, but insufficient biological information for this analysis."""


def read_fasta(path, aligned=False):
    records, name, parts = {}, None, []

    def save():
        if name is None:
            return
        sequence = "".join(parts).upper()
        if not sequence or set(sequence) - (DNA | ({"-"} if aligned else set())):
            raise ValueError(f"Empty sequence or unsupported DNA character: {path}: {name}")
        if name in records or not SAFE_ID.fullmatch(name):
            raise ValueError(f"Duplicate or unsafe FASTA identifier: {path}: {name}")
        records[name] = sequence

    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA header: {path}")
                name, parts = fields[0], []
            elif name is None:
                raise ValueError(f"Sequence before FASTA header: {path}")
            else:
                parts.append(line)
    save()
    return records


def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_table(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def digest_records(records):
    h = hashlib.sha256()
    for name, sequence in sorted(records.items()):
        h.update((name + "\n" + sequence + "\n").encode())
    return h.hexdigest()


def is_true(value):
    if str(value).lower() not in ("true", "false", "1", "0"):
        raise ValueError(f"Invalid boolean in well summary: {value!r}")
    return str(value).lower() in ("true", "1")


def select_consensuses(summary, accepted, review, batch):
    """Join final statuses to the two non-overlapping delivered FASTA files."""
    with Path(summary).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "batch", "well", "status", "diagnostic_available", "diagnostic_kind"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Well summary lacks required final-decision columns")
        rows = list(reader)
    samples = [row["sample"] for row in rows]
    if len(set(samples)) != len(samples):
        raise ValueError("Duplicate samples in well summary")
    ac, rv = read_fasta(accepted), read_fasta(review)
    expected_ac = {r["sample"] for r in rows if r["status"] == "CONSENSUS_PASS"}
    expected_rv = {r["sample"] for r in rows if r["status"] == "REVIEW" and is_true(r["diagnostic_available"])}
    if set(ac) != expected_ac or set(rv) != expected_rv:
        raise ValueError("Delivered FASTA IDs do not match final PASS/available REVIEW summary rows")
    chosen, metadata, excluded = {}, [], []
    batch_rows = [r for r in rows if r["batch"] == batch]
    if not batch_rows:
        raise ValueError(f"Batch absent from final well summary: {batch}")
    for row in batch_rows:
        sample, status = row["sample"], row["status"]
        if sample != f"{batch}__{row['well']}" or not SAFE_ID.fullmatch(sample):
            raise ValueError(f"Sample/batch/well mismatch or unsafe identifier: {sample}")
        if status == "CONSENSUS_PASS":
            kind, sequence = "PASS", ac[sample]
        elif sample in rv:
            kinds = {"DIAGNOSTIC_POLISHED": "REVIEW_POLISHED", "DIAGNOSTIC_UNPOLISHED": "REVIEW_UNPOLISHED"}
            if row["diagnostic_kind"] not in kinds:
                raise ValueError(f"Unknown diagnostic kind for {sample}")
            kind, sequence = kinds[row["diagnostic_kind"]], rv[sample]
        elif status in ("REVIEW", "REJECT", "NO_READS"):
            excluded.append(dict(sample=sample, well=row["well"], status=status,
                                 reason="No available diagnostic consensus" if status == "REVIEW" else "Status excluded",
                                 qc_reasons=row.get("reasons", "")))
            continue
        else:
            raise ValueError(f"Unexpected final well status: {sample}: {status}")
        tip = f"{kind}__{sample}"
        chosen[tip] = sequence
        metadata.append(dict(tip_id=tip, sample=sample, batch=batch, well=row["well"],
                             status=status, kind=kind, construct=row.get("dominant_construct", ""),
                             length=len(sequence), masked_bases=sequence.count("N"),
                             unambiguous_pct=100 * sum(b in "ACGT" for b in sequence) / len(sequence),
                             qc_reasons=row.get("reasons", "")))
    return chosen, metadata, excluded


def validate_alignment(records, aligned):
    if set(records) != set(aligned) or len({len(s) for s in aligned.values()}) != 1:
        raise ValueError("MAFFT output has missing/extra IDs or unequal aligned lengths")
    for name, sequence in records.items():
        if aligned[name].replace("-", "") != sequence:
            raise ValueError(f"MAFFT changed sequence content or orientation: {name}")


def crop_references(refs, aligned, targets, batch, min_overlap):
    """Project the common supported amplicon span through a reference-only MSA."""
    validate_alignment(refs, aligned)
    if targets.get("batch") != batch or targets.get("reference_digest") != digest_records(refs):
        raise ValueError("Target inference batch/reference digest mismatch")
    if set(targets.get("references", {})) != set(refs):
        raise ValueError("Target inference does not describe the complete reference panel")
    donors = {}
    for name, info in targets["references"].items():
        if info.get("status") != "INFERRED" or not info.get("consensus_eligible"):
            continue
        start, end = info.get("start"), info.get("end")
        if (type(start) is not int or type(end) is not int or
                not 0 <= start < end <= len(refs[name])):
            raise ValueError(f"Invalid inferred interval for {name}")
        columns = [i for i, b in enumerate(aligned[name]) if b != "-"]
        donors[name] = (columns[start], columns[end - 1] + 1, start, end)
    if not donors:
        raise NotAnalyzable("No consensus-eligible reference has supported inferred amplicon boundaries")
    left, right = max(x[0] for x in donors.values()), min(x[1] for x in donors.values())
    if left >= right:
        raise NotAnalyzable("Inferred amplicon spans do not overlap in the reference alignment")
    overlap = {name: len(aligned[name][left:right].replace("-", "")) / (end - start)
               for name, (_, _, start, end) in donors.items()}
    if min(overlap.values()) < min_overlap:
        raise NotAnalyzable("Common homologous crop covers too little of at least one inferred target")
    cropped, rows = {}, []
    for name, sequence in aligned.items():
        start = len(sequence[:left].replace("-", ""))
        end = start + len(sequence[left:right].replace("-", ""))
        if start == end:
            raise NotAnalyzable(f"Reference {name} has no bases in the common target span")
        cropped[f"REFERENCE__{name}"] = refs[name][start:end]
        info = targets["references"][name]
        rows.append(dict(reference=name, start=start, end=end, length=end - start,
                         method="common span projected through reference-only MAFFT alignment",
                         own_target_status=info.get("status", ""),
                         own_inferred_start=info.get("start"), own_inferred_end=info.get("end"),
                         crop_donors=";".join(sorted(donors)), msa_start=left, msa_end=right,
                         donor_retained_fraction=overlap.get(name)))
    return cropped, rows


def alignment_metrics(aligned):
    """Retain invariant sites; remove only globally undetermined N/gap columns."""
    sequences = list(aligned.values())
    width = len(sequences[0])
    retained, variable, informative = [], 0, 0
    for i in range(width):
        column = [s[i] for s in sequences]
        if any(b not in "N-" for b in column):
            retained.append(i)
        counts = Counter(b for b in column if b in "ACGT")
        variable += len(counts) > 1
        informative += sum(n >= 2 for n in counts.values()) >= 2
    trimmed = {name: "".join(seq[i] for i in retained) for name, seq in aligned.items()}
    groups = defaultdict(list)
    for name, seq in trimmed.items():
        groups[seq.replace("-", "N")].append(name)
    rows = [{"tip_id": name, "aligned_length": width,
             "unambiguous_pct": 100 * sum(b in "ACGT" for b in seq) / width,
             "ambiguous_pct": 100 * sum(b not in "ACGT-" for b in seq) / width,
             "gap_pct": 100 * seq.count("-") / width} for name, seq in aligned.items()]
    return trimmed, rows, dict(tips=len(aligned), alignment_columns=width,
                              analyzed_columns=len(retained), removed_undetermined_columns=width - len(retained),
                              variable_sites=variable, parsimony_informative_sites=informative,
                              unique_patterns=len(groups),
                              identical_tip_groups=[names for names in groups.values() if len(names) > 1])


def run_command(command, directory, log_name, commands, output=None):
    command = [str(x) for x in command]
    commands.append(command)
    with (directory / log_name).open("w") as log:
        log.write("COMMAND: " + json.dumps(command) + "\n")
        log.flush()
        if output:
            with Path(output).open("w") as handle:
                subprocess.run(command, cwd=directory, stdout=handle, stderr=log, check=True)
        else:
            subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)


def version(executable, directory, name, commands):
    run_command([executable, "--version"], directory, name + "_version.log", commands)
    lines = (directory / (name + "_version.log")).read_text().splitlines()[1:]
    return next((line.strip() for line in lines if re.search(r"\d+\.\d+", line)), "Version output unrecognized; see log")


def align(records, directory, stem, mafft, threads, commands):
    source, destination = directory / f"{stem}.fasta", directory / f"{stem}_aligned.fasta"
    write_fasta(source, records)
    if len(records) == 1:
        write_fasta(destination, records)
    else:
        run_command([mafft, "--auto", "--thread", threads, "--threadit", "0", "--inputorder", source],
                    directory, stem + "_mafft.log", commands, destination)
    result = read_fasta(destination, aligned=True)
    validate_alignment(records, result)
    return result


def inspect_and_render_tree(path, metadata, svg_path):
    from Bio import Phylo
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    tree = Phylo.read(path, "newick")
    tips = tree.get_terminals()
    expected = {row["tip_id"] for row in metadata}
    if len(tips) != len(expected) or {tip.name for tip in tips} != expected:
        raise ValueError("RAxML-NG tree lost, renamed or duplicated input tips")
    if any(c.branch_length is not None and c.branch_length < 0 for c in tree.find_clades()):
        raise ValueError("Negative branch length in inferred tree")
    # Arbitrary display orientation only. No biological rooting or midpoint rooting.
    tree.rooted = False
    tree.ladderize()
    colors = {"PASS": "#168477", "REFERENCE": "#286ca8", "REVIEW_POLISHED": "#bf7c0a", "REVIEW_UNPOLISHED": "#bd4b11"}
    by_tip = {row["tip_id"]: row for row in metadata}
    labels = {tip: (f"REFERENCE | {by_tip[tip]['construct']}" if by_tip[tip]["kind"] == "REFERENCE"
                   else f"{by_tip[tip]['kind']} | {by_tip[tip]['well']} | {by_tip[tip]['construct']}") for tip in expected}
    fig, ax = plt.subplots(figsize=(14, max(4.5, len(tips) * .25 + 1.8)))
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
               label_func=lambda clade: labels.get(clade.name),
               label_colors={labels[name]: colors[by_tip[name]["kind"]] for name in expected},
               branch_labels=lambda clade: f"{clade.confidence:g}" if clade.confidence is not None else None)
    ax.set_title(f"{metadata[0]['batch']} — exploratory unrooted tree\nNode labels: Felsenstein bootstrap support (%)", fontsize=11)
    ax.set_xlabel("Substitutions per site (arbitrary display orientation)")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.legend(handles=[Line2D([], [], color=color, label=kind.replace("_", " ")) for kind, color in colors.items()],
              loc="lower left", bbox_to_anchor=(0, -.08), ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    support = [c.confidence for c in tree.get_nonterminals() if c.confidence is not None]
    return dict(supported_splits=len(support), median_bootstrap_support=statistics.median(support) if support else None,
                splits_below_70=sum(value < 70 for value in support))


def check_tree_tips(path, expected):
    from Bio import Phylo
    count = 0
    for tree in Phylo.parse(path, "newick"):
        tips = [clade.name for clade in tree.get_terminals()]
        if len(tips) != len(expected) or set(tips) != expected:
            raise ValueError(f"Tree lost, renamed or duplicated input tips: {path}")
        count += 1
    return count


def run(args):
    if not SAFE_ID.fullmatch(args.batch):
        raise ValueError("Unsafe batch identifier")
    if min(args.threads, args.bootstraps, args.seed) < 1 or not 0 < args.min_target_overlap <= 1:
        raise ValueError("Threads, bootstraps and seed must be positive; target overlap must be in (0,1]")
    # A nucleotide model string, not a partition/model file. Subprocess uses argv, never a shell.
    if not re.fullmatch(r"[A-Za-z0-9+_.{}=/,-]+", args.model):
        raise ValueError("Invalid nucleotide model string")
    # Resolve relative executable paths before subprocess changes its working directory.
    for name in ("mafft", "raxml_ng"):
        executable = getattr(args, name)
        if "/" in executable:
            setattr(args, name, str(Path(executable).resolve()))
    out = Path(args.outdir).resolve()
    if out.exists():
        raise FileExistsError(f"Output directory already exists: {out}; use a fresh directory")
    out.mkdir(parents=True)
    inputs = {name: Path(getattr(args, name)).resolve() for name in ("references", "targets", "summary", "accepted", "review")}
    report = dict(version=VERSION, batch=args.batch, status="PENDING", reason="", limitations=LIMITATIONS,
                  settings=dict(model=args.model, bootstraps=args.bootstraps, seed=args.seed, threads=args.threads,
                                min_target_overlap=args.min_target_overlap, bootstrap_metric="FBP"),
                  input_sha256={name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in inputs.items()},
                  commands=[], tool_versions={}, counts={}, alignment={}, tree={})
    sequences, metadata, excluded = select_consensuses(inputs["summary"], inputs["accepted"], inputs["review"], args.batch)
    report["counts"] = {kind: sum(r["kind"] == kind for r in metadata) for kind in ("PASS", "REVIEW_POLISHED", "REVIEW_UNPOLISHED")}
    report["counts"].update(references=0, excluded_wells=len(excluded),
                             review_without_sequence=sum(r["status"] == "REVIEW" for r in excluded))
    write_table(out / "excluded_wells.tsv", excluded, ["sample", "well", "status", "reason", "qc_reasons"])
    crops, alignment_rows = [], []
    try:
        refs = read_fasta(inputs["references"])
        if not refs:
            raise ValueError("Reference panel is empty")
        targets = json.loads(inputs["targets"].read_text())
        report["tool_versions"]["mafft"] = version(args.mafft, out, "mafft", report["commands"])
        reference_msa = align(refs, out, "references", args.mafft, args.threads, report["commands"])
        cropped, crops = crop_references(refs, reference_msa, targets, args.batch, args.min_target_overlap)
        write_fasta(out / "cropped_references.fasta", cropped)
        report["counts"]["references"] = len(cropped)
        for tip, sequence in cropped.items():
            metadata.append(dict(tip_id=tip, sample="", batch=args.batch, well="", status="REFERENCE", kind="REFERENCE",
                                 construct=tip.removeprefix("REFERENCE__"), length=len(sequence), masked_bases=sequence.count("N"),
                                 unambiguous_pct=100 * sum(b in "ACGT" for b in sequence) / len(sequence), qc_reasons=""))
        sequences.update(cropped)
        msa = align(sequences, out, "sequences", args.mafft, args.threads, report["commands"])
        tree_msa, alignment_rows, report["alignment"] = alignment_metrics(msa)
        write_fasta(out / "tree_input.fasta", tree_msa)
        if any(not any(b in "ACGT" for b in seq) for seq in tree_msa.values()):
            raise NotAnalyzable("At least one requested tip has no unambiguous aligned bases; no tips were silently dropped")
        if len(tree_msa) < 4 or report["alignment"]["unique_patterns"] < 4:
            raise NotAnalyzable("Fewer than four distinguishable sequences; alignment retained, no tree estimated")
        if not report["alignment"]["parsimony_informative_sites"]:
            raise NotAnalyzable("No parsimony-informative nucleotide sites; alignment retained, no supported tree estimated")
        report["tool_versions"]["raxml_ng"] = version(args.raxml_ng, out, "raxml_ng", report["commands"])
        prefix = out / "raxml" / "tree"
        prefix.parent.mkdir()
        run_command([args.raxml_ng, "--all", "--msa", out / "tree_input.fasta", "--model", args.model,
                     "--bs-trees", args.bootstraps, "--bs-metric", "fbp", "--seed", args.seed,
                     "--threads", f"auto{{{args.threads}}}", "--prefix", prefix],
                    out, "raxml_ng.log", report["commands"])
        # RAxML's optional reduced MSA is deliberately NOT used: identical tips must remain.
        for suffix, destination in (("bestTree", "best_tree.nwk"), ("support", "supported_tree.nwk"),
                                    ("bootstraps", "bootstrap_trees.nwk"), ("bestModel", "fitted_model.txt")):
            source = Path(str(prefix) + ".raxml." + suffix)
            if suffix == "support" and not source.exists():
                source = Path(str(prefix) + ".raxml.supportFBP")
            shutil.copyfile(source, out / destination)
        expected = set(sequences)
        if check_tree_tips(out / "best_tree.nwk", expected) != 1:
            raise ValueError("Expected exactly one best tree")
        n_bootstrap = check_tree_tips(out / "bootstrap_trees.nwk", expected)
        if n_bootstrap != args.bootstraps:
            raise ValueError(f"Expected {args.bootstraps} bootstrap trees; got {n_bootstrap}")
        report["tree"] = inspect_and_render_tree(out / "supported_tree.nwk", metadata, out / "tree.svg")
        report["tree"]["bootstrap_trees"] = n_bootstrap
        report["tree"]["fitted_model"] = (out / "fitted_model.txt").read_text().strip()
        report.update(status="COMPLETE", reason="All requested available consensus and reference tips retained")
    except NotAnalyzable as error:
        report.update(status="SKIPPED", reason=str(error))
    # Missing inputs, malformed provenance, tool failures, and lost tips are hard failures.
    # They must not be presented as a biologically justified skip.
    write_table(out / "sequence_metadata.tsv", metadata, ["tip_id", "sample", "batch", "well", "status", "kind", "construct",
                                                         "length", "masked_bases", "unambiguous_pct", "qc_reasons"])
    write_table(out / "reference_crops.tsv", crops, ["reference", "start", "end", "length", "method", "own_target_status",
                                                   "own_inferred_start", "own_inferred_end", "crop_donors", "msa_start", "msa_end", "donor_retained_fraction"])
    write_table(out / "alignment_metrics.tsv", alignment_rows, ["tip_id", "aligned_length", "unambiguous_pct", "ambiguous_pct", "gap_pct"])
    (out / "phylogeny.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"{args.batch}: {report['status']}: {report['reason']}")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("batch", "references", "targets", "summary", "accepted", "review", "outdir"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--raxml-ng", default="raxml-ng")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--model", default="GTR+G")
    parser.add_argument("--bootstraps", type=int, default=200)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--min-target-overlap", type=float, default=.85)
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
