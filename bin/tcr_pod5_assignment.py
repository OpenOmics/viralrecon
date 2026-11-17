"""Comparison-only informative-position assignment on an aligned reference panel.

Read bases are projected from their primary SAM alignment into a reference-only MSA.
No read is realigned, trimmed, changed or promoted into production consensus evidence.
"""
from collections import Counter
from dataclasses import asdict
import csv
import hashlib
import itertools
import json
from pathlib import Path

from tcr_pod5_core import (DNA, VERSION, Settings, alignment_identity, assign, mean_q,
                           observations, read_fasta)
from tcr_pod5_targets import reference_digest


class InformativePanel:
    def __init__(self, alignment, refs):
        if set(alignment) != set(refs):
            raise ValueError("Aligned reference IDs must exactly match the full competitive panel")
        if len(set(map(len, alignment.values()))) != 1:
            raise ValueError("Reference alignment rows must have equal length")
        for name, seq in alignment.items():
            if seq.replace("-", "") != refs[name]:
                raise ValueError(f"Aligned reference does not reproduce the original sequence: {name}")
        self.alignment = alignment
        self.refs = refs
        self.columns, self.excluded_variable_columns = [], []
        width = len(next(iter(alignment.values())))
        for col in range(width):
            bases = {seq[col] for seq in alignment.values()}
            if len(bases) > 1:
                if bases <= DNA:
                    self.columns.append(col)
                else:
                    self.excluded_variable_columns.append(col)
        self.positions = {}
        for name, seq in alignment.items():
            self.positions[name] = [col for col, base in enumerate(seq) if base != "-"]

    def score(self, records, cfg):
        primary = [r for r in records if not r.flag & (256 | 2048)]
        if len(primary) != 1:
            raise ValueError("Informative scoring requires one primary alignment per read")
        r = primary[0]
        result = dict(assignment=None, best_reference=None, informative_positions=0,
                      deleted_informative_positions=0, best_identity=None, margin=None, identities={})
        if r.flag & 4:
            result["assignment"] = "unmapped"
        elif any(x.flag & 2048 for x in records) or "SA" in r.tags or "N" in r.cigar:
            result["assignment"] = "split_alignment"
        elif mean_q(r.qual) < cfg.min_read_q:
            result["assignment"] = "low_read_quality"
        elif "H" in r.cigar or sum(n for n, op in r.ops if op in "MI=X") / max(1, len(r.seq)) < cfg.min_query_coverage:
            result["assignment"] = "partial_read"
        if result["assignment"]:
            return result
        bases, _ = observations(r, self.refs[r.ref], cfg)
        mapped = {self.positions[r.ref][p]: base for p, base in bases.items()
                  if 0 <= p < len(self.positions[r.ref])}
        observed = {c: mapped[c] for c in self.columns if mapped.get(c) in DNA}
        result["informative_positions"] = len(observed)
        result["deleted_informative_positions"] = sum(mapped.get(c) == "-" for c in self.columns)
        # A deleted diagnostic position supplies no nucleotide identity evidence.
        # All candidates use the SAME covered, canonical columns and denominator.
        identities = {name: sum(seq[c] == base for c, base in observed.items()) / len(observed)
                      for name, seq in self.alignment.items()} if observed else {}
        ranked = sorted(identities, key=lambda name: (-identities[name], name))
        if len(ranked) >= 2:
            best, second = ranked[:2]
            result.update(best_reference=best, best_identity=identities[best],
                          margin=identities[best] - identities[second], identities=identities)
        if len(observed) < cfg.informative_min_positions or len(ranked) < 2:
            result["assignment"] = "insufficient_informative_positions"
        elif result["best_identity"] < cfg.informative_min_identity:
            result["assignment"] = "low_informative_identity"
        elif result["margin"] < cfg.informative_min_margin:
            result["assignment"] = "ambiguous_informative_reference"
        else:
            result["assignment"] = result["best_reference"]
        return result


def evaluate_assignment(a, alignments):
    cfg = Settings.from_json(a.settings)
    refs = read_fasta(a.references)
    panel = InformativePanel(read_fasta(a.reference_alignment, allow_gaps=True), refs)
    manifest = json.loads(Path(a.manifest).read_text())
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=False)
    seen, wells, summaries, transitions = set(), set(), {}, Counter()
    fields = ["read_id", "well", "primary_reference", "legacy_assignment", "deletion_aware_assignment",
              "informative_assignment", "informative_positions", "deleted_informative_positions",
              "best_reference", "best_identity", "margin", "identities", "raw_identity",
              "deletion_aware_identity", "exempt_deleted_bases", "internal_deletions"]
    with open(out / "read_assignment_comparison.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for item in sorted(manifest["alignments"], key=lambda item: item["well"]):
            well = item["well"]
            if well in wells:
                raise ValueError(f"Repeated well in assignment audit: {well}")
            wells.add(well)
            counts = {method: Counter() for method in ("legacy", "deletion_aware", "informative")}
            changed, informative_changed = 0, 0
            for name, group in itertools.groupby(alignments(item["bam"], a.samtools), key=lambda r: r.name):
                if name in seen:
                    raise ValueError("Repeated read UUID; use query-name sorted BAMs with disjoint well IDs")
                seen.add(name)
                records = list(group)
                legacy, _ = assign(records, refs, cfg, require_target_coverage=False, deletion_aware=False)
                current, _ = assign(records, refs, cfg, require_target_coverage=False)
                score = panel.score(records, cfg)
                primary = next(r for r in records if not r.flag & (256 | 2048))
                metrics = alignment_identity(primary, cfg) if not primary.flag & 4 else {}
                informative = score.pop("assignment")
                for method, status in (("legacy", legacy), ("deletion_aware", current), ("informative", informative)):
                    counts[method][status] += 1
                changed += legacy != current
                informative_changed += current != informative
                transitions[(legacy, current, informative)] += 1
                row = dict(read_id=name, well=well, primary_reference=primary.ref,
                           legacy_assignment=legacy, deletion_aware_assignment=current,
                           informative_assignment=informative, **score, **metrics)
                for field in ("identities", "internal_deletions"):
                    row[field] = json.dumps(row.get(field, {} if field == "identities" else []), sort_keys=True)
                writer.writerow(row)
            summaries[well] = dict(counts=counts, deletion_changed_reads=changed,
                                   informative_changed_reads=informative_changed)
    report = dict(software_version=VERSION, batch=manifest.get("batch"), comparison_only=True,
                  reference_digest=reference_digest(refs),
                  alignment_sha256=hashlib.sha256(Path(a.reference_alignment).read_bytes()).hexdigest(),
                  informative_columns_0based=panel.columns,
                  excluded_variable_columns_0based=panel.excluded_variable_columns,
                  settings=asdict(cfg), reads=len(seen), wells=summaries,
                  transitions=[dict(legacy=l, deletion_aware=d, informative=i, reads=n)
                               for (l, d, i), n in sorted(transitions.items())],
                  interpretation=("Comparison only; production assignments and mixture gates are unchanged by this audit. "
                    "Informative scoring uses quality-filtered bases projected through the primary alignment. "
                    "MAPQ/AS and whole-read identity gates are intentionally not applied to this alternative. "
                    "Alignment errors/repeats can misplace bases; inspect changed reads before adopting it. "
                    "Gap/N-containing reference columns and read deletions are not nucleotide identity evidence."))
    (out / "assignment_evaluation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
