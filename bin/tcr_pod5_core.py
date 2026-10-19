#!/usr/bin/env python3
"""Pure-Python validation, SAM evidence and conservative TCR well screening.

This is a research screening heuristic, not a validated single-template assay.
No statistical confidence level or detection limit is implied by the thresholds.
"""
from __future__ import annotations
import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, replace
from itertools import combinations
from pathlib import Path
from typing import Iterable

VERSION = "0.4.0"
DNA = set("ACGT")
WELL = re.compile(r"[A-H](?:[1-9]|1[0-2])$")
SAFE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*$")


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta(path, allow_gaps=False) -> dict[str, str]:
    seqs, name, parts = {}, None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in seqs:
                        raise ValueError(f"Duplicate FASTA ID: {name}")
                    seqs[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name is None:
                raise ValueError("Sequence encountered before FASTA header")
            else:
                parts.append(line)
    if name is not None:
        if name in seqs:
            raise ValueError(f"Duplicate FASTA ID: {name}")
        seqs[name] = "".join(parts).upper()
    alphabet = set("ACGTN-" if allow_gaps else "ACGTN")
    if not seqs or any(not s or set(s) - alphabet for s in seqs.values()):
        raise ValueError("FASTA contains empty records or unsupported nucleotide symbols")
    return seqs


def write_fasta(path, seqs: dict[str, str]):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def edit_distance(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[-1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def full_primer_match(primer: str, window: str):
    """Best alignment of the FULL primer to a substring of a read-end window.

    Query-global, target-local edit distance. Returns (errors, start, end),
    with 0-based, half-open coordinates. Not used to assign well barcodes.
    """
    prev = [(0, j) for j in range(len(window) + 1)]
    for i, x in enumerate(primer, 1):
        cur = [(i, 0)]
        for j, y in enumerate(window, 1):
            choices = [(prev[j][0] + 1, prev[j][1]),
                       (cur[-1][0] + 1, cur[-1][1]),
                       (prev[j - 1][0] + (x != y), prev[j - 1][1])]
            cur.append(min(choices, key=lambda t: (t[0], -t[1])))
        prev = cur
    end = min(range(1, len(prev)), key=lambda j: (prev[j][0], j))
    return prev[end][0], prev[end][1], end


def parse_primers(path, expected_pairs=48):
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"Shorthand", "Direction", "Sequence", "I-start", "I-end", "Pair_with"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Primer TSV requires columns: {sorted(required)}")
        rows = [{k: (v or "").strip() for k, v in r.items() if k is not None}
                for r in reader if any((v or "").strip() for v in r.values() if isinstance(v, str))]
    by_name = {}
    for r in rows:
        name = r["Shorthand"]
        if name in by_name:
            raise ValueError(f"Duplicate primer shorthand: {name}")
        seq = r["Sequence"].upper()
        if not seq or set(seq) - DNA:
            raise ValueError(f"Non-ACGT primer: {name}")
        r["Sequence"] = seq
        try:
            start, end = int(r["I-start"]) - 1, int(r["I-end"])
        except ValueError as e:
            raise ValueError(f"Invalid index coordinates for {name}") from e
        if not 0 <= start < end <= len(seq):
            raise ValueError(f"Index coordinates outside primer: {name}")
        if r["Direction"] not in {"F", "R"}:
            raise ValueError(f"Invalid direction for {name}")
        r.update(prefix=seq[:start], index=seq[start:end], body=seq[end:])
        if not r["body"]:
            raise ValueError(f"Missing post-index primer body: {name}")
        by_name[name] = r
    pairs, used = [], set()
    for name, f in by_name.items():
        if f["Direction"] != "F":
            continue
        if not name.startswith("f") or not WELL.fullmatch(name[1:]):
            raise ValueError(f"Expected fA1-style shorthand, received {name}")
        well = name[1:]
        mate = f["Pair_with"]
        if mate != "r" + well or mate not in by_name:
            raise ValueError(f"Missing or inconsistent reverse mate for {name}: {mate}")
        r = by_name[mate]
        if r["Direction"] != "R" or (r["Pair_with"] and r["Pair_with"] != name):
            raise ValueError(f"Invalid reciprocal pairing: {name}, {mate}")
        if mate in used:
            raise ValueError(f"Reused reverse primer: {mate}")
        used.add(mate)
        pairs.append({"well": well, "forward": f, "reverse": r})
    unused = [n for n, r in by_name.items() if r["Direction"] == "R" and n not in used]
    if unused:
        raise ValueError(f"Unpaired reverse primers: {unused}")
    if len(pairs) != expected_pairs:
        raise ValueError(f"Expected {expected_pairs} complete primer pairs; found {len(pairs)}")
    pairs.sort(key=lambda p: (p["well"][0], int(p["well"][1:])))
    for side in ("forward", "reverse"):
        for field in ("prefix", "body"):
            if len({p[side][field] for p in pairs}) != 1:
                raise ValueError(f"Dorado arrangement requires shared {side} {field}")
        indexes = [p[side]["index"] for p in pairs]
        if len(set(indexes)) != len(indexes):
            raise ValueError(f"Duplicate {side} indexes")
    if len({len(p[s]["index"]) for p in pairs for s in ("forward", "reverse")}) != 1:
        raise ValueError("Custom barcode sequences must all have the same length")
    return pairs


def build_scheme(primers, refs, outdir, expected_pairs=48, barcode_errors=1, expected_constructs=None):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=False)
    (out / "primers").mkdir()
    pairs = parse_primers(primers, expected_pairs)
    if barcode_errors < 0 or barcode_errors >= len(pairs[0]["forward"]["index"]):
        raise ValueError("Invalid barcode error limit")
    distances = {}
    for side in ("forward", "reverse"):
        ds = [edit_distance(a[side]["index"], b[side]["index"])
              for a, b in combinations(pairs, 2)]
        distances[side] = min(ds) if ds else None
    # A conservative start. No assumption that the real assay tolerates these values.
    f, r = pairs[0]["forward"], pairs[0]["reverse"]
    toml = f'''[arrangement]
name = "TCR_WELLS"
kit = "TCR"
mask1_front = "{f['prefix']}"
mask1_rear = "{f['body']}"
mask2_front = "{r['prefix']}"
mask2_rear = "{r['body']}"
barcode1_pattern = "FBC%02i"
barcode2_pattern = "RBC%02i"
first_index = 1
last_index = {len(pairs)}

[scoring]
max_barcode_penalty = {barcode_errors}
min_barcode_penalty_dist = 2
min_separation_only_dist = 100
flank_left_pad = 0
flank_right_pad = 0
front_barcode_window = 150
rear_barcode_window = 150
barcode_end_proximity = 100
min_flank_score = 0.8
midstrand_flank_score = 0.95
'''
    (out / "barcodes.toml").write_text(toml)
    bcseqs = {}
    for i, p in enumerate(pairs, 1):
        p["barcode_number"] = i
        p["barcode_label"] = f"TCR_WELLS_barcode{i:02d}"
        bcseqs[f"FBC{i:02d}"] = p["forward"]["index"]
        bcseqs[f"RBC{i:02d}"] = p["reverse"]["index"]
        write_fasta(out / "primers" / (p["well"] + ".fasta"), {
            f"{p['well']}_front\tet:Z:primer\tsk:Z:any": p["forward"]["Sequence"],
            f"{p['well']}_rear\tet:Z:primer\tsk:Z:any": p["reverse"]["Sequence"],
        })
    write_fasta(out / "barcodes.fasta", bcseqs)
    records = read_fasta(refs)
    if len(records) > 50:
        raise ValueError("This implementation supports at most 50 competing TCR constructs")
    if any(not SAFE.fullmatch(name) for name in records):
        raise ValueError("Reference FASTA IDs must contain only letters, digits, underscores, periods and hyphens")
    if len(set(records.values())) != len(records):
        raise ValueError("Reference panel contains identical sequences; assignment is not identifiable")
    expected = sorted(records if expected_constructs is None else expected_constructs)
    if not expected or len(set(expected)) != len(expected) or set(expected) - set(records):
        raise ValueError("expected_constructs must be a nonempty, unique subset of reference FASTA IDs")
    # Preserve ALL full reference sequences, including contamination sentinels.
    write_fasta(out / "reference.fasta", records)
    (out / "scheme.json").write_text(json.dumps({
        "pairs": pairs, "index_coordinates": "1-based inclusive",
        "index_min_edit_distances": distances,
        "expected_constructs": expected, "reference_lengths": {k: len(v) for k, v in records.items()},
        "reference_mode": "full_panel", "version": VERSION,
        "threshold_status": "exploratory, not assay-validated"
    }, indent=2) + "\n")
    with open(out / "well_map.tsv", "w") as ftable:
        ftable.write("well\tbarcode\tforward\treverse\n")
        for p in pairs:
            ftable.write(f"{p['well']}\t{p['barcode_label']}\t{p['forward']['index']}\t{p['reverse']['index']}\n")
    return pairs


@dataclass
class Settings:
    min_reads: int = 50
    polish_min_reads: int = 20
    polish_min_callable_fraction: float = 0.80
    min_baseq: int = 15
    min_read_q: float = 10.0
    min_identity: float = 0.90
    deletion_min_length: int = 50
    deletion_min_anchor: int = 20
    informative_min_identity: float = 0.75
    informative_min_margin: float = 0.15
    informative_min_positions: int = 50
    diagnostic_min_reads: int = 10
    diagnostic_min_depth: int = 10
    diagnostic_min_strand_reads: int = 2
    diagnostic_min_callable_fraction: float = 0.10
    min_query_coverage: float = 0.85
    min_target_coverage: float = 0.85
    min_mapq: int = 20
    score_gap: int = 20
    minor_fraction: float = 0.10
    minor_reads: int = 5
    minor_each_strand: int = 2
    alert_fraction: float = 0.03
    alert_reads: int = 3
    max_unresolved_fraction: float = 0.10
    max_trim_fail_fraction: float = 0.10
    min_strand_reads: int = 5
    min_depth: int = 20
    consensus_fraction: float = 0.80
    min_callable_fraction: float = 0.98
    phase_spacing: int = 10
    primer_errors: int = 3
    primer_window: int = 150
    max_overtrim: int = 5
    homopolymer_length: int = 5
    target_min_reads: int = 50
    target_min_wells: int = 2
    target_min_reads_per_well: int = 5
    target_min_each_strand: int = 5
    target_min_alignment_length: int = 3000
    target_min_query_coverage: float = 0.95
    target_boundary_tolerance: int = 75
    target_min_cluster_fraction: float = 0.80

    @classmethod
    def from_json(cls, path):
        raw = json.loads(Path(path).read_text())
        unknown = set(raw) - set(asdict(cls()))
        if unknown:
            raise ValueError(f"Unknown QC options: {sorted(unknown)}")
        obj = cls(**raw)
        for key, value in asdict(obj).items():
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid numeric QC option {key}={value}")
            if isinstance(getattr(cls(), key), int) and not isinstance(value, int):
                raise ValueError(f"QC option {key} must be an integer")
        for key in ("min_identity", "min_query_coverage", "min_target_coverage", "minor_fraction",
                    "alert_fraction", "max_unresolved_fraction", "max_trim_fail_fraction",
                    "consensus_fraction", "min_callable_fraction", "polish_min_callable_fraction",
                    "informative_min_identity", "informative_min_margin", "diagnostic_min_callable_fraction",
                    "target_min_query_coverage", "target_min_cluster_fraction"):
            if not 0 < getattr(obj, key) <= 1:
                raise ValueError(f"QC option {key} must be in (0,1]")
        if not obj.alert_fraction < obj.minor_fraction <= 0.5:
            raise ValueError("Require alert_fraction < minor_fraction <= 0.5")
        if obj.consensus_fraction <= 0.5 or obj.min_reads < 1 or obj.polish_min_reads < 1 or obj.min_depth < 1:
            raise ValueError("Invalid consensus or support thresholds")
        for key in ("target_min_reads", "target_min_wells", "target_min_reads_per_well",
                    "target_min_each_strand", "target_min_alignment_length", "target_boundary_tolerance",
                    "deletion_min_length", "deletion_min_anchor", "informative_min_positions",
                    "diagnostic_min_reads", "diagnostic_min_depth", "diagnostic_min_strand_reads"):
            if getattr(obj, key) < 1:
                raise ValueError(f"QC option {key} must be positive")
        if obj.target_min_cluster_fraction <= 0.5:
            raise ValueError("target_min_cluster_fraction must be greater than 0.5")
        return obj


@dataclass
class Sam:
    name: str
    flag: int
    ref: str
    start: int
    mapq: int
    cigar: str
    seq: str
    qual: str
    tags: dict[str, str]
    line: str

    @classmethod
    def parse(cls, line):
        a = line.rstrip("\n").split("\t")
        if len(a) < 11:
            raise ValueError("Malformed SAM alignment")
        tags = {}
        for field in a[11:]:
            p = field.split(":", 2)
            if len(p) == 3:
                tags[p[0]] = p[2]
        if a[9] != "*" and a[10] != "*" and len(a[9]) != len(a[10]):
            raise ValueError(f"SAM sequence/quality length mismatch: {a[0]}")
        return cls(a[0], int(a[1]), a[2], int(a[3]) - 1, int(a[4]), a[5], a[9], a[10], tags, line)

    @property
    def strand(self):
        return "-" if self.flag & 16 else "+"

    @property
    def ops(self):
        if self.cigar == "*":
            return []
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", self.cigar)]
        if "".join(f"{n}{op}" for n, op in ops) != self.cigar:
            raise ValueError(f"Malformed CIGAR: {self.cigar}")
        return ops

    @property
    def end(self):
        return self.start + sum(n for n, op in self.ops if op in "MDN=X")


def mean_q(qual):
    if qual == "*" or not qual:
        return 0.0
    return -10 * math.log10(sum(10 ** (-(ord(c) - 33) / 10) for c in qual) / len(qual))


def alignment_identity(r: Sam, cfg: Settings):
    """Keep long, internally anchored deletions separate from base-error identity.

    Only D operations with contiguous aligned sequence on BOTH sides qualify.
    Short/unanchored deletions, insertions and mismatches still count as errors.
    This does not establish that the deletion is real; evidence() tests its support.
    """
    ops, deleted, pos = r.ops, [], r.start
    span = sum(n for n, op in ops if op in "MID=X")
    for i, (n, op) in enumerate(ops):
        if op == "D" and n >= cfg.deletion_min_length:
            anchors = []
            for step in (-1, 1):
                j, length = i + step, 0
                while 0 <= j < len(ops) and ops[j][1] in "M=X":
                    length += ops[j][0]
                    j += step
                anchors.append(length)
            if min(anchors) >= cfg.deletion_min_anchor:
                deleted.append((pos, n))
        if op in "MDN=X":
            pos += n
    exempt = sum(n for _, n in deleted)
    # Missing NM is never silently interpreted as a perfect alignment.
    errors = int(r.tags.get("NM", span))
    if not 0 <= errors <= span or errors < sum(n for n, op in ops if op in "ID"):
        raise ValueError(f"Inconsistent NM/CIGAR for read {r.name}")
    matches = span - errors
    return dict(raw_identity=matches / max(1, span),
                deletion_aware_identity=matches / max(1, span - exempt),
                exempt_deleted_bases=exempt, internal_deletions=deleted)


def assign(records: list[Sam], refs: dict[str, str], cfg: Settings, require_target_coverage=True,
           deletion_aware=True):
    primary = [r for r in records if not r.flag & (256 | 2048)]
    if len(primary) != 1:
        raise ValueError(f"Expected one primary alignment per read: {records[0].name}")
    r = primary[0]
    if r.flag & 4:
        return "unmapped", None
    if r.ref not in refs:
        raise ValueError(f"Unexpected reference in SAM: {r.ref}")
    if any(x.flag & 2048 for x in records) or "SA" in r.tags:
        return "split_alignment", None
    if mean_q(r.qual) < cfg.min_read_q:
        return "low_read_quality", None
    if "H" in r.cigar:
        return "partial_read", None
    if "N" in r.cigar:
        return "split_alignment", None
    aligned = sum(n for n, op in r.ops if op in "MI=X")
    metrics = alignment_identity(r, cfg)
    identity = metrics["deletion_aware_identity" if deletion_aware else "raw_identity"]
    if identity < cfg.min_identity:
        return "low_identity", None
    if aligned / max(1, len(r.seq)) < cfg.min_query_coverage:
        return "partial_read", None
    if require_target_coverage and (r.end - r.start) / len(refs[r.ref]) < cfg.min_target_coverage:
        return "partial_target", None
    alt = [int(x.tags.get("AS", -10**9)) for x in records
           if x.ref != r.ref and not x.flag & (4 | 2048)]
    gap = int(r.tags.get("AS", -10**9)) - max(alt) if alt else 10**9
    if r.mapq == 255 or r.mapq < cfg.min_mapq or gap < cfg.score_gap:
        return "ambiguous_reference", None
    return r.ref, r


def normalize_event(kind, pos, seq, ref):
    # I: position is the preceding reference base. D: first deleted base.
    if kind == "I":
        while pos >= 0 and ref[pos] == seq[-1]:
            seq = seq[-1] + seq[:-1]
            pos -= 1
    else:
        while pos > 0 and ref[pos - 1] == seq[-1]:
            seq = seq[-1] + seq[:-1]
            pos -= 1
    return kind, pos, seq


def observations(r: Sam, ref: str, cfg: Settings):
    q, p, bases, events, deletions = 0, r.start, {}, {}, []
    if r.qual == "*":
        return bases, events
    qs = [ord(x) - 33 for x in r.qual]
    for n, op in r.ops:
        if op in "M=X":
            for k in range(n):
                if 0 <= p + k < len(ref) and qs[q + k] >= cfg.min_baseq and r.seq[q + k] in DNA:
                    bases[p + k] = r.seq[q + k]
            p, q = p + n, q + n
        elif op == "I":
            iq = qs[max(0, q - 1):min(len(qs), q + n + 1)]
            if q > 0 and q + n < len(qs) and iq and min(iq) >= cfg.min_baseq:
                ev = normalize_event("I", p - 1, r.seq[q:q + n], ref)
                events[ev] = True
            q += n
        elif op == "D":
            if q > 0 and q < len(qs) and min(qs[q - 1], qs[q]) >= cfg.min_baseq:
                ev = normalize_event("D", p, ref[p:p + n], ref)
                events[ev] = True
                deletions.append((p, n, ev[1]))
                for k in range(n):
                    bases[p + k] = "-"
            p += n
        elif op == "N":
            p += n
        elif op == "S":
            q += n
        elif op in "HP":
            pass
    # Project equivalent repeat deletions onto the same coordinates as their
    # normalized event. Otherwise the draft deletes one position but masks another.
    for original, length, normalized in deletions:
        if normalized < original:
            shifted = {p + length: bases.get(p) for p in range(normalized, original)}
            for p, base in shifted.items():
                if base in DNA:
                    bases[p] = base
                else:
                    bases.pop(p, None)
            for p in range(normalized, normalized + length):
                bases[p] = "-"
    return bases, events


def low_complexity(ref, pos, n):
    s = ref[max(0, pos - n):min(len(ref), pos + n + 1)]
    return any(base * n in s for base in DNA)


def evidence(reads: list[Sam], ref: str, cfg: Settings, interval=None, mask_variation=False):
    """Evaluate bases in a fixed target interval; report original-reference coordinates."""
    start, end = interval if interval is not None else (0, len(ref))
    if not 0 <= start < end <= len(ref):
        raise ValueError("Invalid evidence interval")
    obs, evs, strands = {}, {}, {}
    pile = defaultdict(dict)
    all_events = set()
    for r in reads:
        b, e = observations(r, ref, cfg)
        if interval is not None:
            b = {pos: base for pos, base in b.items() if start <= pos < end}
            e = {event: yes for event, yes in e.items()
                 if (start <= event[1] and event[1] + (len(event[2]) if event[0] == "D" else 1) < end)}
        obs[r.name], evs[r.name], strands[r.name] = b, e, r.strand
        all_events.update(e)
        for pos, base in b.items():
            pile[pos][r.name] = base
    sites = []
    for pos, calls in pile.items():
        calls = {k: v for k, v in calls.items() if v != "-"}
        counts = Counter(calls.values())
        if len(counts) >= 2:
            ranked = counts.most_common()
            minor = ranked[1][1]
            if minor >= cfg.alert_reads and minor / sum(counts.values()) >= cfg.alert_fraction:
                sites.append(dict(kind="SNP", pos=pos, alleles=dict(counts), calls=calls,
                                  homopolymer=low_complexity(ref, pos, cfg.homopolymer_length)))
    event_support = []
    for event in sorted(all_events):
        kind, pos, seq = event
        calls = {}
        # Require flanking base coverage for reference/no-event observations.
        check = [pos, pos + 1] if kind == "I" else list(range(pos - 1, pos + len(seq) + 1))
        for rid, bases in obs.items():
            if event in evs[rid]:
                calls[rid] = "ALT"
            elif all(bases.get(p) in DNA for p in check):
                nearby = any(abs(e[1] - pos) <= max(len(seq), len(e[2])) + 1 for e in evs[rid])
                if not nearby:
                    calls[rid] = "REF"
        counts = Counter(calls.values())
        event_support.append((event, calls))
        if len(counts) == 2:
            minor = min(counts.values())
            if minor >= cfg.alert_reads and minor / sum(counts.values()) >= cfg.alert_fraction:
                sites.append(dict(kind=kind, pos=pos, sequence=seq, alleles=dict(counts), calls=calls,
                                  homopolymer=(not (kind == "D" and len(seq) >= cfg.deletion_min_length)
                                               and low_complexity(ref, pos, cfg.homopolymer_length))))
    # Annotate every sensitive candidate site with the stronger, strand-balanced
    # evidence used for consensus blocking. Candidate detection remains intentionally
    # sensitive (alert_reads / alert_fraction); acceptance uses the stronger rule.
    for site in sites:
        counts = Counter(site["alleles"])
        ranked = counts.most_common()
        minor_allele, minor_count = ranked[1]
        informative = sum(counts.values())
        minor_ids = [rid for rid, allele in site["calls"].items() if allele == minor_allele]
        minor_strands = Counter(strands[rid] for rid in minor_ids)
        strong_minor = (minor_count >= cfg.minor_reads
                        and minor_count / max(1, informative) >= cfg.minor_fraction
                        and min(minor_strands.get("+", 0), minor_strands.get("-", 0)) >= cfg.minor_each_strand)
        site.update(minor_allele=minor_allele, minor_count=minor_count,
                    minor_fraction=minor_count / max(1, informative),
                    minor_strand_counts=dict(minor_strands), strong_minor=strong_minor,
                    blocking_variation=strong_minor and not site["homopolymer"])

    links = []
    for a, b in combinations(sites, 2):
        if a["homopolymer"] or b["homopolymer"] or abs(a["pos"] - b["pos"]) < cfg.phase_spacing:
            continue
        groups = defaultdict(list)
        for rid in a["calls"].keys() & b["calls"].keys():
            groups[(a["calls"][rid], b["calls"][rid])].append(rid)
        n = sum(map(len, groups.values()))
        strong = []
        for hap, ids in groups.items():
            st = Counter(strands[rid] for rid in ids)
            if (len(ids) >= cfg.minor_reads and len(ids) / max(1, n) >= cfg.minor_fraction
                    and min(st.get("+", 0), st.get("-", 0)) >= cfg.minor_each_strand):
                strong.append((hap, ids))
        for (h1, ids1), (h2, ids2) in combinations(strong, 2):
            if h1[0] != h2[0] and h1[1] != h2[1]:
                links.append(dict(site1=a["pos"], site2=b["pos"],
                                  haplotype1=h1, haplotype2=h2,
                                  reads1=len(ids1), reads2=len(ids2), informative_reads=n))
    # A conservative draft, based on high-quality reads only, never reference fill.
    draft, callable_count = [], 0
    forced_mask = set()
    if mask_variation:
        for site in sites:
            if site.get("strong_minor"):
                length = len(site["sequence"]) if site["kind"] == "D" else (2 if site["kind"] == "I" else 1)
                forced_mask.update(range(site["pos"], site["pos"] + length))
    fixed = []
    for event, calls in event_support:
        counts = Counter(calls.values())
        kind, pos, sequence = event
        affected = set(range(pos, pos + (len(sequence) if kind == "D" else 2)))
        if (not affected & forced_mask and sum(counts.values()) >= cfg.min_depth
                and counts["ALT"] / sum(counts.values()) >= cfg.consensus_fraction):
            fixed.append(event)
    deleted, insertions, conflicts = set(), {}, False
    for kind, pos, seq in fixed:
        if kind == "D":
            positions = set(range(pos, pos + len(seq)))
            if deleted & positions:
                conflicts = True
            deleted.update(positions)
        else:
            if pos in insertions and insertions[pos] != seq:
                conflicts = True
            insertions[pos] = seq
    if set(insertions) & deleted:
        conflicts = True
    if -1 in insertions:
        draft.append(insertions[-1])
    for p in range(start, end):
        counts = Counter(pile.get(p, {}).values())
        depth = sum(counts.values())
        base, n = counts.most_common(1)[0] if counts else ("N", 0)
        called = p not in forced_mask and depth >= cfg.min_depth and n / max(1, depth) >= cfg.consensus_fraction
        if p in deleted:
            if called and base == "-":
                callable_count += 1
            else:
                conflicts = True
        else:
            if called and base in DNA:
                draft.append(base)
                callable_count += 1
            else:
                draft.append("N")
        if p in insertions:
            draft.append(insertions[p])
    return dict(sites=sites, links=links, pile=pile, events=event_support,
                draft="".join(draft), callable_fraction=callable_count / (end - start),
                conflicting_fixed_indels=conflicts, strands=strands,
                masked_variation_positions=sorted(forced_mask))


def screen(assignments, assigned_reads, refs, cfg: Settings, trim_stats=None,
           targets=None, expected_constructs=None):
    """Return (report, evidence, selected reads). Never claim physical template count."""
    counts = Counter(assignments)
    construct_counts = {k: counts.get(k, 0) for k in refs}
    ranked = sorted(construct_counts, key=lambda k: (-construct_counts[k], k))
    dominant = ranked[0] if construct_counts[ranked[0]] else None
    n_assigned = sum(construct_counts.values())
    total = sum(counts.values())
    reasons, warnings, strong_mixture = [], [], []
    if total < cfg.min_reads or construct_counts.get(dominant, 0) < cfg.min_reads:
        reasons.append("LOW_SUPPORT")
    for ref in ranked[1:]:
        n = construct_counts[ref]
        fraction = n / max(1, n_assigned)
        if n >= cfg.alert_reads and fraction >= cfg.alert_fraction:
            st = Counter(r.strand for r in assigned_reads if r.ref == ref)
            if (n >= cfg.minor_reads and fraction >= cfg.minor_fraction and
                    min(st.get("+", 0), st.get("-", 0)) >= cfg.minor_each_strand):
                strong_mixture.append(ref)
            else:
                reasons.append("REVIEW_SECOND_CONSTRUCT")
    if strong_mixture:
        reasons.append("REJECT_MULTIPLE_CONSTRUCTS")

    # Shared vector sequence can make reads non-discriminating between otherwise
    # valid references. Keep those reads out of consensus evidence, but do not count
    # them as failed/unresolved reads for the unresolved-fraction gate.
    ambiguous_reference_reads = counts.get("ambiguous_reference", 0)
    unresolved_for_qc = total - n_assigned - ambiguous_reference_reads
    if ambiguous_reference_reads:
        warnings.append("WARN_AMBIGUOUS_REFERENCE")
    if unresolved_for_qc / max(1, total) > cfg.max_unresolved_fraction:
        reasons.append("REVIEW_UNRESOLVED_READS")
    if counts.get("split_alignment", 0) >= cfg.alert_reads:
        reasons.append("REVIEW_SPLIT_ALIGNMENTS")
    if trim_stats:
        n, failed = trim_stats["input_reads"], trim_stats["failed_reads"]
        if failed / max(1, n) > cfg.max_trim_fail_fraction:
            warnings.append("WARN_PRIMER_TRIMMING")

    selected = [r for r in assigned_reads if r.ref == dominant]
    target = (targets.get(dominant) if targets is not None else (0, len(refs[dominant]))) if dominant else None
    target_partial_reads = 0
    target_boundary_outlier_reads = 0
    if expected_constructs is not None and dominant and dominant not in expected_constructs:
        reasons.append("REVIEW_UNEXPECTED_CONSTRUCT")
    if target is None:
        reasons.append("REVIEW_TARGET_NOT_INFERRED")
    else:
        start, end = target
        if targets is not None:
            target_boundary_outlier_reads = sum(abs(r.start - start) > cfg.target_boundary_tolerance or
                                                abs(r.end - end) > cfg.target_boundary_tolerance for r in selected)
            if target_boundary_outlier_reads / max(1, len(selected)) > cfg.max_unresolved_fraction:
                reasons.append("REVIEW_TARGET_BOUNDARIES")
        complete = [r for r in selected
                    if max(0, min(end, r.end) - max(start, r.start)) / (end - start) >= cfg.min_target_coverage]
        target_partial_reads = len(selected) - len(complete)
        if target_partial_reads / max(1, len(selected)) > cfg.max_unresolved_fraction:
            reasons.append("REVIEW_PARTIAL_TARGET_READS")
        selected = complete
        if len(selected) < cfg.min_reads:
            reasons.append("LOW_TARGET_SUPPORT")
    strands = Counter(r.strand for r in selected)
    if min(strands.get("+", 0), strands.get("-", 0)) < cfg.min_strand_reads:
        reasons.append("REVIEW_STRAND_SUPPORT")

    if target is None:
        ev = dict(sites=[], links=[], pile={}, events=[], draft="", callable_fraction=None,
                  conflicting_fixed_indels=False, strands=dict(strands))
    else:
        ev = evidence(selected, refs[dominant], cfg, target)
    strong_sites = [site for site in ev["sites"] if site.get("strong_minor")]
    blocking_sites = [site for site in strong_sites if site.get("blocking_variation")]
    if ev["links"]:
        reasons.append("REJECT_LINKED_HAPLOTYPES")
    elif blocking_sites:
        reasons.append("REVIEW_WITHIN_CONSTRUCT_VARIATION")
    elif ev["sites"]:
        warnings.append("WARN_WITHIN_CONSTRUCT_VARIATION")
    if any(site.get("strong_minor") and site.get("homopolymer") for site in ev["sites"]):
        warnings.append("WARN_HOMOPOLYMER_VARIATION")
    if ev["callable_fraction"] is not None and ev["callable_fraction"] < cfg.min_callable_fraction:
        reasons.append("REVIEW_INCOMPLETE_COVERAGE")
    if ev["conflicting_fixed_indels"]:
        reasons.append("REVIEW_CONFLICTING_INDELS")

    if any(x.startswith("REJECT") for x in reasons):
        status = "REJECT"
    elif reasons:
        status = "REVIEW"
    else:
        status = "PASS_SCREEN"

    # Polishing is a reconstruction step, not an acceptance decision. Some REVIEW
    # wells may enter polish when there is enough complete, strand-balanced dominant
    # evidence and no hard/structural blocker. Final acceptance remains stricter.
    polish_blockers = []
    persistent = {
        "REVIEW_SECOND_CONSTRUCT", "REVIEW_SPLIT_ALIGNMENTS", "REVIEW_UNRESOLVED_READS",
        "REVIEW_UNEXPECTED_CONSTRUCT", "REVIEW_TARGET_NOT_INFERRED", "REVIEW_TARGET_BOUNDARIES",
        "REVIEW_PARTIAL_TARGET_READS", "REVIEW_CONFLICTING_INDELS",
    }
    for reason in reasons:
        if reason.startswith("REJECT") or reason in persistent:
            polish_blockers.append(reason)
    if len(selected) < cfg.polish_min_reads:
        polish_blockers.append("POLISH_LOW_SUPPORT")
    if min(strands.get("+", 0), strands.get("-", 0)) < cfg.min_strand_reads:
        polish_blockers.append("POLISH_STRAND_SUPPORT")
    if ev["callable_fraction"] is None or ev["callable_fraction"] < cfg.polish_min_callable_fraction:
        polish_blockers.append("POLISH_INCOMPLETE_COVERAGE")
    polish_eligible = not polish_blockers

    report = dict(status=status, reasons=sorted(set(reasons)), warnings=sorted(set(warnings)),
                  polish_eligible=polish_eligible, polish_blockers=sorted(set(polish_blockers)),
                  dominant_construct=dominant, total_reads=total, assigned_reads=n_assigned,
                  selected_reads=len(selected),
                  construct_counts=construct_counts,
                  unresolved_counts={k: v for k, v in counts.items() if k not in refs},
                  ambiguous_reference_reads=ambiguous_reference_reads,
                  unresolved_for_qc=unresolved_for_qc,
                  strand_counts=dict(strands), callable_fraction=ev["callable_fraction"],
                  candidate_sites=len(ev["sites"]), strong_candidate_sites=len(strong_sites),
                  blocking_variation_sites=len(blocking_sites), linked_haplotypes=ev["links"],
                  target_interval=list(target) if target is not None else None,
                  target_partial_reads=target_partial_reads,
                  target_boundary_outlier_reads=target_boundary_outlier_reads,
                  evidence_coordinates="original reference, 0-based; draft spans target_interval",
                  expected_constructs=sorted(expected_constructs) if expected_constructs is not None else None,
                  settings=asdict(cfg), software_version=VERSION,
                  interpretation=("Polish eligibility permits reconstruction of some REVIEW wells; "
                                  "CONSENSUS_PASS remains a separate post-polish acceptance decision."))
    return report, ev, selected


def diagnostic_consensus(report, reads, refs, cfg: Settings):
    """Masked, unpolished draft only. Never changes the well's acceptance decision."""
    blockers = [r for r in report["reasons"] if r.startswith("REJECT")]
    for reason in ("REVIEW_TARGET_NOT_INFERRED", "REVIEW_UNEXPECTED_CONSTRUCT", "REVIEW_CONFLICTING_INDELS"):
        if reason in report["reasons"]:
            blockers.append(reason)
    if report["status"] == "REJECT" and not blockers:
        blockers.append("REJECT")
    if not report.get("dominant_construct") or not report.get("target_interval"):
        blockers.append("DIAGNOSTIC_NO_TARGET")
    if len(reads) < cfg.diagnostic_min_reads:
        blockers.append("DIAGNOSTIC_LOW_SUPPORT")
    strands = Counter(r.strand for r in reads)
    if min(strands.get("+", 0), strands.get("-", 0)) < cfg.diagnostic_min_strand_reads:
        blockers.append("DIAGNOSTIC_STRAND_SUPPORT")
    ev = None
    if not blockers:
        ev = evidence(reads, refs[report["dominant_construct"]],
                      replace(cfg, min_depth=cfg.diagnostic_min_depth),
                      report["target_interval"], mask_variation=True)
        if ev["conflicting_fixed_indels"]:
            blockers.append("DIAGNOSTIC_CONFLICTING_INDELS")
        if ev["callable_fraction"] < cfg.diagnostic_min_callable_fraction or not set(ev["draft"]) & DNA:
            blockers.append("DIAGNOSTIC_INCOMPLETE_COVERAGE")
    detail = dict(eligible=not blockers, blockers=sorted(set(blockers)),
                  kind="DIAGNOSTIC_UNPOLISHED", accepted=False, selected_reads=len(reads),
                  strand_counts=dict(strands), min_depth=cfg.diagnostic_min_depth,
                  min_baseq=cfg.min_baseq, base_support=cfg.consensus_fraction,
                  callable_fraction=ev["callable_fraction"] if ev else None,
                  source_reasons=report["reasons"], source_warnings=report["warnings"],
                  interpretation="Dominant-construct subset only; masked draft, not a validated single-template call")
    return detail, ev["draft"] if ev and not blockers else None
