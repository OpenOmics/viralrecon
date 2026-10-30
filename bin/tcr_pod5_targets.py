#!/usr/bin/env python3
"""Batch-level, well-balanced empirical target inference (no user BED required).

Intervals are estimates in ORIGINAL reference coordinates, not independently
validated PCR boundaries. Never substitute a full reference for a failed inference.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import hashlib
import math
from statistics import median

from tcr_pod5_core import Sam, Settings, VERSION


@dataclass(frozen=True)
class Endpoint:
    read_id: str
    well: str
    reference: str
    start: int
    end: int
    strand: str


def reference_digest(refs):
    """Sequence-based digest, independent of FASTA wrapping and record order."""
    h = hashlib.sha256()
    for name, seq in sorted(refs.items()):
        h.update((name + "\n" + seq + "\n").encode())
    return h.hexdigest()


def target_candidate(read: Sam, cfg: Settings):
    """Additional inference-only filters, AFTER competitive assignment.

    Short reads still count as construct evidence in well QC, but cannot define
    boundaries. Target coverage cannot be used here: it has not been inferred yet.
    """
    aligned = sum(n for n, op in read.ops if op in "MI=X")
    if aligned < cfg.target_min_alignment_length:
        return "short_alignment"
    if aligned / max(1, len(read.seq)) < cfg.target_min_query_coverage:
        return "clipped_alignment"
    return None


def quantile(values, fraction):
    values = sorted(values)
    index = (len(values) - 1) * fraction
    lo, hi = math.floor(index), math.ceil(index)
    return values[lo] + (values[hi] - values[lo]) * (index - lo)


def endpoint_summary(reads):
    return {
        axis: {"p05": quantile(values, .05), "median": median(values), "p95": quantile(values, .95)}
        for axis, values in (("start", [r.start for r in reads]), ("end", [r.end for r in reads]))
    } if reads else {}


def infer_one(reference, length, observations, cfg):
    by_well = defaultdict(list)
    for r in observations:
        if not (0 <= r.start < r.end <= length):
            raise ValueError(f"Alignment outside reference {reference}: {r}")
        by_well[r.well].append(r)
    eligible = {w: rs for w, rs in sorted(by_well.items()) if len(rs) >= cfg.target_min_reads_per_well}
    report = dict(reference=reference, reference_length=length, status="INSUFFICIENT_SUPPORT",
                  start=None, end=None, candidate_reads=len(observations), candidate_wells=len(by_well),
                  eligible_reads=sum(map(len, eligible.values())), eligible_wells=len(eligible),
                  supporting_reads=0, supporting_wells=0, balanced_cluster_fraction=0.0,
                  strand_counts={}, all_endpoint_summary=endpoint_summary(observations),
                  per_well_candidates={w: len(rs) for w, rs in sorted(by_well.items())})
    if report["eligible_reads"] < cfg.target_min_reads or len(eligible) < cfg.target_min_wells:
        report["reason"] = "Insufficient long, high-quality reads from independent wells"
        return report, set()

    tolerance = cfg.target_boundary_tolerance
    candidates = sorted(set((median(r.start for r in rs), median(r.end for r in rs))
                            for rs in eligible.values()))

    def cluster(center):
        # Each well contributes equal total weight regardless of sequencing depth.
        selected = {w: [r for r in rs if abs(r.start - center[0]) <= tolerance
                       and abs(r.end - center[1]) <= tolerance] for w, rs in eligible.items()}
        score = sum(len(selected[w]) / len(eligible[w]) for w in eligible) / len(eligible)
        return selected, score

    ranked = []
    for center in candidates:
        selected, score = cluster(center)
        n_wells = sum(len(rs) >= cfg.target_min_reads_per_well for rs in selected.values())
        n_reads = sum(map(len, selected.values()))
        ranked.append((score, n_wells, n_reads, center, selected))
    # Deterministic ties; a comparable second population fails the fraction gate.
    score, n_wells, n_reads, center, selected = sorted(
        ranked, key=lambda x: (-x[0], -x[1], -x[2], x[3]))[0]
    supported = {w: rs for w, rs in selected.items() if len(rs) >= cfg.target_min_reads_per_well}
    if not supported:
        report.update(status="UNSTABLE_BOUNDARIES", reason="No supported joint endpoint cluster")
        return report, set()
    center = (median(median(r.start for r in rs) for rs in supported.values()),
              median(median(r.end for r in rs) for rs in supported.values()))
    selected, score = cluster(center)
    supported = {w: rs for w, rs in selected.items() if len(rs) >= cfg.target_min_reads_per_well}
    support = [r for rs in supported.values() for r in rs]
    strands = Counter(r.strand for r in support)
    report.update(supporting_reads=len(support), supporting_wells=len(supported),
                  balanced_cluster_fraction=score, strand_counts=dict(strands),
                  cluster_endpoint_summary=endpoint_summary(support),
                  per_well_support={w: len(rs) for w, rs in sorted(supported.items())})
    if score < cfg.target_min_cluster_fraction:
        report.update(status="UNSTABLE_BOUNDARIES", reason="No dominant joint start/end population across wells")
    elif len(support) < cfg.target_min_reads or len(supported) < cfg.target_min_wells:
        report["reason"] = "Dominant boundary population has insufficient read or well support"
    elif min(strands.get("+", 0), strands.get("-", 0)) < cfg.target_min_each_strand:
        report["reason"] = "Dominant boundary population lacks support on both alignment strands"
    else:
        # Round outwards when medians fall between integer reference boundaries.
        start, end = math.floor(center[0]), math.ceil(center[1])
        if end - start < cfg.target_min_alignment_length:
            report.update(status="UNSTABLE_BOUNDARIES", reason="Inferred reference span is shorter than the configured minimum")
        else:
            report.update(status="INFERRED", reason="Well-balanced, strand-supported joint endpoint cluster",
                          start=start, end=end, length=end - start)
            return report, {r.read_id for r in support}
    return report, set()


def infer_targets(refs, observations, cfg, batch, expected_constructs=None):
    expected = set(refs) if expected_constructs is None else set(expected_constructs)
    if not expected or expected - set(refs):
        raise ValueError("Expected constructs are not a nonempty subset of the reference panel")
    grouped, seen = defaultdict(list), set()
    for row in observations:
        if row.reference not in refs:
            raise ValueError(f"Unknown reference {row.reference}")
        if row.read_id in seen:
            raise ValueError(f"Repeated read UUID in target inference: {row.read_id}")
        seen.add(row.read_id)
        grouped[row.reference].append(row)
    results, used = {}, set()
    for name in sorted(refs):
        results[name], selected = infer_one(name, len(refs[name]), grouped[name], cfg)
        results[name]["consensus_eligible"] = name in expected
        used.update(selected)
    return dict(batch=batch, version=VERSION, reference_digest=reference_digest(refs),
                coordinate_system="original reference; 0-based start-inclusive, end-exclusive",
                method="well-balanced joint endpoint cluster; median of per-well medians",
                expected_constructs=sorted(expected), settings=asdict(cfg), references=results,
                limitations="Empirical boundaries are not independent proof of a complete PCR product; systematic truncation or shared deletions can be missed."), used


def load_intervals(report, refs, batch):
    """Check provenance before applying a batch's intervals to any well."""
    if report.get("batch") != batch or report.get("reference_digest") != reference_digest(refs):
        raise ValueError("Target inference batch/reference does not match this well")
    if set(report.get("references", {})) != set(refs):
        raise ValueError("Target inference must report every reference, including unsupported sentinels")
    expected = report.get("expected_constructs", [])
    if not expected or set(expected) - set(refs):
        raise ValueError("Invalid expected-construct annotation in target report")
    intervals = {}
    for ref, row in report["references"].items():
        if row["status"] != "INFERRED":
            continue
        start, end = row["start"], row["end"]
        if not isinstance(start, int) or not isinstance(end, int) or not 0 <= start < end <= len(refs[ref]):
            raise ValueError(f"Invalid inferred interval for {ref}")
        intervals[ref] = (start, end)
    return intervals
