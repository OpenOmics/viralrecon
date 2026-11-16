"""Synthetic regression tests; these are not sequencing results."""
import argparse
import copy
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
import tcr_pod5_phylogeny as phy


class PhylogenyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.batch = "batch_one"
        self.refs = {"construct_A": "TTACGTACGTAA", "construct_B": "TTACGTTCGTAA",
                     "sentinel_C": "TTTCGTACGTAA", "sentinel_D": "TTTCGTTCGTAA"}
        phy.write_fasta(self.root / "refs.fa", self.refs)
        self.targets = dict(batch=self.batch, reference_digest=phy.digest_records(self.refs), references={
            name: dict(status="INFERRED" if i < 2 else "INSUFFICIENT_SUPPORT", consensus_eligible=i < 2,
                       start=2 if i < 2 else None, end=10 if i < 2 else None)
            for i, name in enumerate(self.refs)})
        (self.root / "targets.json").write_text(json.dumps(self.targets))
        self.rows = []
        for well, status, kind, available in (("A1", "CONSENSUS_PASS", "", False),
                                             ("A2", "REVIEW", "DIAGNOSTIC_POLISHED", True),
                                             ("A3", "REVIEW", "DIAGNOSTIC_UNPOLISHED", True),
                                             ("A4", "REVIEW", "", False), ("A5", "REJECT", "", False),
                                             ("A6", "NO_READS", "", False)):
            self.rows.append(dict(sample=f"{self.batch}__{well}", batch=self.batch, well=well,
                                  status=status, diagnostic_kind=kind, diagnostic_available=available,
                                  dominant_construct="construct_A", reasons="test_reason"))
        # A second batch must never be pulled into this batch's tree.
        self.rows.append(dict(sample="other__B1", batch="other", well="B1", status="CONSENSUS_PASS",
                              diagnostic_kind="", diagnostic_available=False, dominant_construct="construct_B", reasons=""))
        self.save_summary()
        phy.write_fasta(self.root / "accepted.fa", {f"{self.batch}__A1": "ACGTACGT", "other__B1": "ACGTTCGT"})
        phy.write_fasta(self.root / "review.fa", {f"{self.batch}__A2": "ACGTACGT", f"{self.batch}__A3": "NCGTACGT"})
        self.args = argparse.Namespace(batch=self.batch, references=str(self.root / "refs.fa"),
            targets=str(self.root / "targets.json"), summary=str(self.root / "summary.tsv"),
            accepted=str(self.root / "accepted.fa"), review=str(self.root / "review.fa"),
            outdir=str(self.root / "out"), mafft="mafft", raxml_ng="raxml-ng", threads=4,
            bootstraps=2, seed=17, model="GTR+G", min_target_overlap=.85)

    def save_summary(self):
        phy.write_table(self.root / "summary.tsv", self.rows, list(self.rows[0]))

    def select(self):
        return phy.select_consensuses(self.args.summary, self.args.accepted, self.args.review, self.batch)

    def test_final_selection_includes_each_review_once_and_preserves_masks(self):
        sequences, rows, excluded = self.select()
        self.assertEqual(len(sequences), 3)
        self.assertEqual([r["kind"] for r in rows], ["PASS", "REVIEW_POLISHED", "REVIEW_UNPOLISHED"])
        self.assertEqual(sequences["REVIEW_UNPOLISHED__batch_one__A3"], "NCGTACGT")
        self.assertEqual({row["well"] for row in excluded}, {"A4", "A5", "A6"})
        self.assertFalse(any("other" in key for key in sequences))

    def test_missing_promised_sequence_is_error(self):
        phy.write_fasta(self.root / "review.fa", {f"{self.batch}__A2": "ACGTACGT"})
        with self.assertRaisesRegex(ValueError, "FASTA IDs"):
            self.select()

    def test_rejected_sequence_cannot_slip_into_review_fasta(self):
        with (self.root / "review.fa").open("a") as handle:
            handle.write(">batch_one__A5\nACGTACGT\n")
        with self.assertRaisesRegex(ValueError, "FASTA IDs"):
            self.select()

    def test_duplicate_summary_and_fasta_ids_fail(self):
        self.rows.append(self.rows[0])
        self.save_summary()
        with self.assertRaisesRegex(ValueError, "Duplicate samples"):
            self.select()
        (self.root / "duplicate.fa").write_text(">a\nAC\n>a description\nAC\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            phy.read_fasta(self.root / "duplicate.fa")

    def test_unknown_status_and_diagnostic_kind_fail(self):
        self.rows[1]["diagnostic_kind"] = "UNKNOWN"
        self.save_summary()
        with self.assertRaisesRegex(ValueError, "diagnostic kind"):
            self.select()
        self.rows[1]["diagnostic_kind"] = "DIAGNOSTIC_POLISHED"
        self.rows[3]["status"] = "PASS_SCREEN"
        self.save_summary()
        with self.assertRaisesRegex(ValueError, "Unexpected final"):
            self.select()

    def test_all_panel_references_are_cropped_even_without_own_targets(self):
        cropped, rows = phy.crop_references(self.refs, self.refs, self.targets, self.batch, .85)
        self.assertEqual(len(cropped), 4)
        self.assertEqual(cropped["REFERENCE__sentinel_C"], "TCGTACGT")
        self.assertEqual([(r["start"], r["end"]) for r in rows], [(2, 10)] * 4)
        self.assertEqual(rows[2]["own_target_status"], "INSUFFICIENT_SUPPORT")

    def test_crop_projection_handles_insertions_not_raw_coordinate_reuse(self):
        # X carries an insertion before the target and another inside it.
        aligned = {"A": "--TTAC-GTAA", "X": "CCTTACGGTAA"}
        refs = {name: seq.replace("-", "") for name, seq in aligned.items()}
        targets = dict(batch="b", reference_digest=phy.digest_records(refs), references={
            "A": dict(status="INFERRED", consensus_eligible=True, start=2, end=6),
            "X": dict(status="INSUFFICIENT_SUPPORT", consensus_eligible=False)})
        cropped, rows = phy.crop_references(refs, aligned, targets, "b", .85)
        self.assertEqual(cropped, {"REFERENCE__A": "ACGT", "REFERENCE__X": "ACGGT"})
        self.assertEqual((rows[1]["start"], rows[1]["end"]), (4, 9))

    def test_common_span_guard_and_nonoverlap_skip(self):
        targets = copy.deepcopy(self.targets)
        targets["references"]["construct_B"].update(start=8, end=12)
        with self.assertRaisesRegex(phy.NotAnalyzable, "too little"):
            phy.crop_references(self.refs, self.refs, targets, self.batch, .85)
        targets["references"]["construct_B"].update(start=10, end=12)
        with self.assertRaisesRegex(phy.NotAnalyzable, "do not overlap"):
            phy.crop_references(self.refs, self.refs, targets, self.batch, .85)

    def test_reference_digest_and_invalid_coordinates_fail(self):
        targets = copy.deepcopy(self.targets)
        targets["reference_digest"] = "wrong"
        with self.assertRaisesRegex(ValueError, "digest mismatch"):
            phy.crop_references(self.refs, self.refs, targets, self.batch, .85)
        targets = copy.deepcopy(self.targets)
        targets["references"]["construct_A"]["start"] = -1
        with self.assertRaisesRegex(ValueError, "Invalid inferred"):
            phy.crop_references(self.refs, self.refs, targets, self.batch, .85)

    def test_alignment_rejects_modified_sequences_or_lost_tips(self):
        with self.assertRaisesRegex(ValueError, "content or orientation"):
            phy.validate_alignment({"A": "ACGT"}, {"A": "TGCA"})
        with self.assertRaisesRegex(ValueError, "missing/extra"):
            phy.validate_alignment({"A": "ACGT", "B": "ACGT"}, {"A": "ACGT"})

    def test_invariant_columns_kept_and_only_undetermined_columns_removed(self):
        msa = {"a": "ACGN-", "b": "ACGN-", "c": "ATGN-", "d": "ATGN-"}
        trimmed, rows, counts = phy.alignment_metrics(msa)
        self.assertEqual(trimmed["a"], "ACG")
        self.assertEqual(counts["analyzed_columns"], 3)
        self.assertEqual(counts["variable_sites"], 1)
        self.assertEqual(counts["parsimony_informative_sites"], 1)
        self.assertEqual(counts["unique_patterns"], 2)
        self.assertEqual(counts["identical_tip_groups"], [["a", "b"], ["c", "d"]])
        self.assertEqual(rows[0]["ambiguous_pct"], 20)
        self.assertEqual(rows[0]["gap_pct"], 20)

    def fake_tools(self):
        """Exercise real subprocess/file boundaries, but never claim real inference."""
        mafft = self.root / "mock_mafft"
        mafft.write_text("#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\n"
                         "print('v7.526 TEST DOUBLE' if '--version' in sys.argv else Path(sys.argv[-1]).read_text(), end='\\n')\n")
        raxml = self.root / "mock_raxml"
        raxml.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
args = sys.argv[1:]
if '--version' in args:
    print('RAxML-NG v1.2.2 TEST DOUBLE')
    sys.exit()
prefix = args[args.index('--prefix') + 1]
names = [line[1:].split()[0] for line in Path(args[args.index('--msa') + 1]).read_text().splitlines() if line.startswith('>')]
tree = '(' + ','.join(name + ':0.01' for name in names) + ');\\n'
for suffix in ('bestTree', 'support'):
    Path(prefix + '.raxml.' + suffix).write_text(tree)
Path(prefix + '.raxml.bestModel').write_text('GTR+G4')
Path(prefix + '.raxml.bootstraps').write_text(tree * int(args[args.index('--bs-trees') + 1]))
''')
        for executable in (mafft, raxml):
            executable.chmod(0o755)
        self.args.mafft, self.args.raxml_ng = str(mafft), str(raxml)

    def test_tool_file_contract_and_no_silent_tip_deduplication(self):
        self.fake_tools()
        def inspect(path, metadata, svg):
            text = Path(path).read_text()
            self.assertTrue(all(row["tip_id"] + ":" in text for row in metadata))
            self.assertEqual(len(metadata), 7)  # three consensus + complete four-reference panel
            Path(svg).write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
            return {"supported_splits": 0}
        with patch.object(phy, "inspect_and_render_tree", side_effect=inspect), \
             patch.object(phy, "check_tree_tips", side_effect=lambda path, expected: Path(path).read_text().count(";")):
            report = phy.run(self.args)
        self.assertEqual(report["status"], "COMPLETE")
        self.assertEqual(report["tree"]["bootstrap_trees"], 2)
        self.assertEqual(report["counts"]["review_without_sequence"], 1)
        command = report["commands"][-1]
        self.assertEqual(command[command.index("--threads") + 1], "auto{4}")
        self.assertNotIn("--force", command)
        self.assertTrue((self.root / "out/phylogeny.json").is_file())
        self.assertEqual(len(phy.read_fasta(self.root / "out/tree_input.fasta", aligned=True)), 7)
        # A second attempt cannot overwrite the evidence from an earlier attempt.
        with self.assertRaises(FileExistsError):
            phy.run(self.args)

    def test_missing_boundary_is_reported_skip_not_full_reference_fallback(self):
        self.fake_tools()
        for info in self.targets["references"].values():
            info["status"] = "INSUFFICIENT_SUPPORT"
        (self.root / "targets.json").write_text(json.dumps(self.targets))
        report = phy.run(self.args)
        self.assertEqual(report["status"], "SKIPPED")
        self.assertFalse((self.root / "out/best_tree.nwk").exists())
        self.assertFalse((self.root / "out/cropped_references.fasta").exists())
        self.assertTrue((self.root / "out/sequence_metadata.tsv").exists())

    def test_external_failure_is_not_converted_to_biological_skip(self):
        self.fake_tools()
        Path(self.args.mafft).write_text("#!/bin/sh\nexit 23\n")
        with self.assertRaises(subprocess.CalledProcessError):
            phy.run(self.args)
        self.assertFalse((self.root / "out/phylogeny.json").exists())

    @unittest.skipUnless(importlib.util.find_spec("Bio"), "Biopython is not installed")
    def test_newick_tip_validation_and_svg_rendering(self):
        sequences, metadata, _ = self.select()
        path = self.root / "tree.nwk"
        names = list(sequences)
        path.write_text(f"(({names[0]}:0.01,{names[1]}:0.02)77:0.1,{names[2]}:0.03);\n")
        result = phy.inspect_and_render_tree(path, metadata, self.root / "tree.svg")
        self.assertEqual(result["median_bootstrap_support"], 77)
        self.assertIn("PASS", (self.root / "tree.svg").read_text())
        self.assertEqual(phy.check_tree_tips(path, set(sequences)), 1)
        with self.assertRaisesRegex(ValueError, "lost, renamed or duplicated"):
            phy.check_tree_tips(path, set(sequences) | {"missing"})

    @unittest.skipUnless(shutil.which("mafft") and shutil.which("raxml-ng") and importlib.util.find_spec("Bio"),
                         "Real MAFFT/RAxML-NG/Biopython integration requires the phylogeny environment")
    def test_real_tools_smoke(self):
        # Deterministic 400 bp fixture with shared informative and private mutations.
        import random
        rng = random.Random(13)
        base = [rng.choice("ACGT") for _ in range(400)]
        refs = {}
        for i, name in enumerate(self.refs):
            seq = base.copy()
            for pos in range(0, 40, 2):
                seq[pos] = "AC"[i // 2]
                seq[pos + 1] = "GT"[i % 2]
            refs[name] = "".join(seq)
        phy.write_fasta(self.root / "refs.fa", refs)
        self.targets["reference_digest"] = phy.digest_records(refs)
        for info in self.targets["references"].values():
            if info["status"] == "INFERRED":
                info.update(start=0, end=400)
        (self.root / "targets.json").write_text(json.dumps(self.targets))
        phy.write_fasta(self.root / "accepted.fa", {"batch_one__A1": refs["construct_A"], "other__B1": refs["construct_B"]})
        phy.write_fasta(self.root / "review.fa", {"batch_one__A2": refs["construct_A"], "batch_one__A3": "N" * 10 + refs["construct_A"][10:]})
        report = phy.run(self.args)
        self.assertEqual(report["status"], "COMPLETE")
        self.assertEqual(report["alignment"]["tips"], 7)


if __name__ == "__main__":
    unittest.main()
