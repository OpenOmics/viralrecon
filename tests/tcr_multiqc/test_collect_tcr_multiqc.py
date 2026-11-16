"""Synthetic regression fixtures; none of these counts are experimental results."""
import argparse
import csv
import json
import re
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "assets/tcr_pod5/multiqc"))
from collect_tcr_multiqc import Report, histogram_points


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


class CollectorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "results" / "tcr_pod5"
        self.args = argparse.Namespace(results_dir=str(self.root.parent), batch=None,
            outdir=str(Path(self.tmp.name) / "report"), audit_dir=[], assignment_evaluation=[],
            no_discover_audits=True, base_config=None, render=False, multiqc="multiqc")
        self.summary = []
        for batch, ns in (("batch_one", [70, 30, 0]), ("batch_two", [20])):
            manifest = [{"batch": batch, "well": f"A{i+1}", "reads": n} for i, n in enumerate(ns)]
            manifest.append({"batch": batch, "well": "unclassified", "reads": 100})
            put(self.root / batch / "demultiplexing/demux/manifest.json", manifest)
            put(self.root / batch / "basecalling/basecall_provenance.json", {"reads": sum(ns) + 100, "dorado_version": "TEST"})
            put(self.root / batch / "target_inference/targets/inferred_targets.json",
                {"batch": batch, "settings": {"min_reads": 50, "min_depth": 20, "min_callable_fraction": .98}, "references": {}})
            for i, n in enumerate(ns):
                well = f"A{i+1}"
                wd = self.root / batch / "wells" / well
                accepted = batch == "batch_one" and i == 0
                status = "CONSENSUS_PASS" if accepted else "REVIEW" if n else "NO_READS"
                selected = n - 15 if n else 0
                row = {"sample": batch + "__" + well, "batch": batch, "well": well, "raw_assigned_reads": n,
                       "status": status, "usable_reads": n - 10 if n else 0, "selected_reads": selected,
                       "callable_fraction": .85 if n else "", "diagnostic_available": status == "REVIEW",
                       "diagnostic_kind": "DIAGNOSTIC_POLISHED" if status == "REVIEW" else "",
                       "diagnostic_callable_fraction": .95 if status == "REVIEW" else ""}
                self.summary.append(row)
                if n:
                    put(wd / "trimming/trim_stats.json", {"input_reads": n, "passed_reads": n-5, "failed_reads": 5, "reasons": {"pass": n-5, "test_failure": 5}})
                    put(wd / "heterogeneity/qc/decision.json", {
                        "status": "REVIEW", "assigned_reads": n-10, "selected_reads": selected, "callable_fraction": .85,
                        "strand_counts": {"+": selected-2, "-": 2}, "reasons": ["test_reason"], "warnings": [],
                        "polish_eligible": True, "polish_blockers": [], "construct_counts": {"ref": n-10},
                        "unresolved_counts": {"low_identity": 5}, "diagnostic": {"eligible": True, "callable_fraction": .4}})
                    put(wd / "consensus/consensus_result/decision.json", {
                        "status": status, "length": 100, "masked_bases": 0 if accepted else 5,
                        "selected_reads": selected-1, "strand_counts": {"+": selected-3, "-": 2}, "reasons": []})
        sf = self.root / "summary/well_summary.tsv"
        sf.parent.mkdir(parents=True)
        with sf.open("w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(self.summary[0]), delimiter="\t")
            wr.writeheader(); wr.writerows(self.summary)
        put(self.root / "summary/run_summary.json", {"wells": 4, "accepted_consensuses": 1, "diagnostic_consensuses": 2, "review_consensuses": 2,
            "statuses": {"CONSENSUS_PASS": 1, "REVIEW": 2, "NO_READS": 1}})
        (self.root / "summary/accepted_consensus.fasta").write_text(">batch_one__A1\n" + "A"*100 + "\n")
        diag = "".join(">" + s + "\n" + "A"*95 + "N"*5 + "\n" for s in ("batch_one__A2", "batch_two__A1"))
        for name in ("diagnostic_consensus.fasta", "review_consensus.fasta"):
            (self.root / "summary" / name).write_text(diag)

    def build(self):
        r = Report(self.args)
        r.load(); r.overview(); r.inputs_and_basecalls(); r.demux_trim_assignment()
        r.target_and_qc(); r.consensus(); r.phylogeny(); r.audits(); r.completeness()
        return r

    def add_phylogeny(self, status="COMPLETE"):
        folder = self.root / "batch_one/phylogeny"
        put(folder / "phylogeny.json", dict(batch="batch_one", version="1.0.0", status=status,
            reason="test reason", limitations="Exploratory; REVIEW remains unaccepted.",
            counts=dict(PASS=1, REVIEW_POLISHED=1, REVIEW_UNPOLISHED=0, references=4,
                        excluded_wells=1, review_without_sequence=0),
            alignment=dict(tips=6, analyzed_columns=100, identical_tip_groups=[["tip1", "tip2"]]),
            settings=dict(model="GTR+G", bootstraps=200, seed=17, threads=4, min_target_overlap=.85),
            tool_versions=dict(mafft="TEST", raxml_ng="TEST"), tree=dict(bootstrap_trees=200) if status == "COMPLETE" else {}))
        for filename, header, row in (
            ("reference_crops.tsv", "reference\tstart\tend\tlength", "ref1\t512\t4136\t3624"),
            ("sequence_metadata.tsv", "tip_id\tkind\tqc_reasons", "tip1\tREVIEW_POLISHED\tREJECT_LINKED_HAPLOTYPES"),
            ("excluded_wells.tsv", "sample\tstatus\treason", "batch_one__A3\tNO_READS\tStatus excluded"),
            ("alignment_metrics.tsv", "tip_id\tunambiguous_pct\tambiguous_pct\tgap_pct", "tip1\t90\t5\t5")):
            (folder / filename).write_text(header + "\n" + row + "\n")
        if status == "COMPLETE":
            (folder / "tree.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>TEST TREE</text></svg>')
            for name in ("best_tree.nwk", "supported_tree.nwk"):
                (folder / name).write_text("(tip1:0.1,tip2:0.1);\n")
            for name in ("tree_input.fasta", "sequences_aligned.fasta"):
                (folder / name).write_text(">tip1\nACGT\n")
        return folder

    def test_phylogeny_is_absent_when_not_requested(self):
        r = self.build()
        self.assertFalse(any(key.startswith(("tcr_17_", "tcr_18_")) for key in r.data))

    def test_phylogeny_tables_plot_downloads_and_original_status(self):
        self.add_phylogeny()
        r = self.build()
        self.assertEqual(r.data["tcr_17_selection"]["data"]["batch_one"]["references"], 4)
        self.assertEqual(r.data["tcr_17_crops"]["data"]["batch_one / ref1"]["start"], 512)
        self.assertEqual(r.data["tcr_17_batch_one_composition"]["data"]["tip1"]["Unambiguous bases"], 90.0)
        page = r.data["tcr_18_batch_one_tree"]["data"]
        self.assertIn("data:image/svg+xml;base64,", page)
        self.assertIn('download="batch_one_supported_tree.nwk"', page)
        self.assertEqual(r.data["tcr_10_final"]["data"]["batch_one / A2"]["status"], "REVIEW")
        r.write()
        self.assertTrue((Path(self.args.outdir) / "plots/batch_one_tree.svg").is_file())

    def test_phylogeny_respects_batch_filter(self):
        self.add_phylogeny()
        self.args.batch = ["batch_two"]
        r = self.build()
        self.assertNotIn("tcr_17_selection", r.data)

    def test_skipped_phylogeny_has_reason_and_no_fake_tree(self):
        self.add_phylogeny(status="SKIPPED")
        r = self.build()
        self.assertEqual(r.data["tcr_18_settings"]["data"]["batch_one"]["status"], "SKIPPED")
        page = r.data["tcr_18_batch_one_tree"]["data"]
        self.assertIn("Tree not estimated: test reason", page)
        self.assertNotIn("data:image", page)

    def test_incomplete_completed_phylogeny_fails(self):
        folder = self.add_phylogeny()
        (folder / "supported_tree.nwk").unlink()
        with self.assertRaisesRegex(ValueError, "Completed phylogeny is missing"):
            self.build()

    def test_display_alias_does_not_rewrite_source_decisions(self):
        path = self.root / "batch_one/wells/A1/consensus/consensus_result/decision.json"
        post = json.loads(path.read_text())
        post["reasons"] = ["REJECT_LINKED_HAPLOTYPES", "REJECT_POSTPOLISH_LINKED_HAPLOTYPES"]
        put(path, post)
        before = path.read_bytes()
        qc_path = self.root / "batch_one/wells/A1/heterogeneity/qc/decision.json"
        qc = json.loads(qc_path.read_text())
        qc["reasons"] = ["REJECT_LINKED_HAPLOTYPES"]
        put(qc_path, qc)
        r = self.build()
        self.assertEqual(path.read_bytes(), before)
        row = r.tables["tcr_09_polished"]["batch_one / A1"]
        self.assertEqual(row["reasons"], "REJECT_LINKED_VARIANTS; REJECT_POSTPOLISH_LINKED_VARIANTS")
        self.assertIn("Consensus outputs and sequences for review", r.data["tcr_10_final"]["parent_name"])
        self.assertEqual(r.data["tcr_08_reasons"]["data"]["batch_one"]["REJECT_LINKED_VARIANTS"], 1)
        self.assertNotIn("REJECT_LINKED_HAPLOTYPES", r.data["tcr_08_reasons"]["headers"])

    def test_full_counts_and_metric_stages(self):
        r = self.build()
        top = r.data["tcr_00_batch"]["data"]
        self.assertEqual(top["batch_one"]["raw_assigned_reads"], 100)
        self.assertEqual(top["batch_one"]["trimmed_reads"], 90)
        self.assertEqual(top["batch_one"]["panel_assigned_reads"], 80)
        self.assertEqual(top["batch_one"]["selected_reads"], 70)
        self.assertEqual(top["batch_one"]["accepted_consensuses"], 1)
        self.assertEqual(top["batch_one"]["diagnostic_consensuses"], 1)
        overview = r.data["tcr_00_well"]["data"]["batch_one / A1"]
        self.assertEqual((overview["status"], overview["selected_reads"], overview["callable_pct"]), ("CONSENSUS_PASS", 55, 85))
        final = r.data["tcr_09_polished"]["data"]["batch_one / A1"]
        self.assertEqual((final["polish_selected_reads"], final["post_callable_pct"]), (54, 100))
        self.assertEqual(r.warnings, [])

    def test_missing_not_zero_and_zero_read_not_missing(self):
        (self.root / "batch_one/wells/A2/trimming/trim_stats.json").unlink()
        r = self.build()
        self.assertNotIn("trimmed_reads", r.data["tcr_00_batch"]["data"]["batch_one"])
        self.assertEqual(r.availability["batch_one / A2"]["trim"], "MISSING")
        self.assertEqual(r.availability["batch_one / A3"]["trim"], "NOT EXPECTED")

    def test_resolved_variation_is_visible_without_reclassifying_reports(self):
        path = self.root / "batch_one/wells/A1/consensus/consensus_result/decision.json"
        post = json.loads(path.read_text())
        post.update(resolved_reasons=["REVIEW_WITHIN_CONSTRUCT_VARIATION"],
                    variation_reassessment={"status": "RESOLVED"})
        put(path, post)
        r = self.build()
        rows = r.data["tcr_09_polished"]["data"]
        self.assertEqual(rows["batch_one / A1"]["resolved_reasons"], "REVIEW_WITHIN_CONSTRUCT_VARIATION")
        self.assertEqual(rows["batch_one / A1"]["variation_reassessment"], "RESOLVED")
        self.assertEqual(rows["batch_one / A1"]["status"], post["status"])
        self.assertEqual(rows["batch_one / A2"]["variation_reassessment"], "")
        self.assertEqual(rows["batch_one / A2"]["status"], "REVIEW")

    def test_split_counts_are_stage_specific_and_missing_values_stay_blank(self):
        wd = self.root / "batch_one/wells/A1"
        for path, count, total in ((wd / "heterogeneity/qc/decision.json", 3, 1000),
                                   (wd / "consensus/consensus_result/decision.json", 4, 200)):
            report = json.loads(path.read_text())
            report["split_alignment_qc"] = {"reads": count, "total_reads": total, "fraction": count / total}
            put(path, report)
        r = self.build()
        for section, count, total, percent in (("tcr_08_qc", 3, 1000, .3),
                                              ("tcr_09_polished", 4, 200, 2.0)):
            rows = r.data[section]["data"]
            self.assertEqual(rows["batch_one / A1"]["split_alignment_reads"], count)
            self.assertEqual(rows["batch_one / A1"]["split_alignment_total_reads"], total)
            self.assertEqual(rows["batch_one / A1"]["split_alignment_pct"], percent)
            self.assertIsNone(rows["batch_one / A2"].get("split_alignment_pct"))
        self.assertEqual(r.data["tcr_09_polished"]["data"]["batch_one / A1"]["status"], "CONSENSUS_PASS")

    def test_batch_selection_and_whole_run_totals(self):
        self.args.batch = ["batch_two"]
        r = self.build()
        self.assertEqual(list(r.data["tcr_00_batch"]["data"]), ["batch_two"])
        self.assertEqual(r.data["tcr_10_totals"]["data"]["Whole run (all batches)"]["wells"], 4)

    def test_duplicate_well_fails(self):
        fp = self.root / "batch_one/demultiplexing/demux/manifest.json"
        d = json.loads(fp.read_text()); d.append(d[0]); put(fp, d)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.build()

    def test_subset_consistency_check(self):
        (self.root / "summary/review_consensus.fasta").write_text(">wrong_sample\nACGT\n")
        r = self.build()
        self.assertTrue(any("identical-sequence subset" in x for x in r.warnings))

    def test_output_and_overwrite_protection(self):
        r = self.build(); r.write()
        out = Path(self.args.outdir)
        self.assertTrue((out / "multiqc_config.yaml").is_file())
        data = json.loads((out / "report_data.json").read_text())
        self.assertEqual(data["custom_content"]["order"][0], "tcr_00")
        self.assertGreater(data["report_section_order"]["tcr_00_batch"]["order"], data["report_section_order"]["tcr_00_well"]["order"])
        with self.assertRaises(FileExistsError):
            r.write()

    def test_histogram_numeric_order_and_tail_accounting(self):
        # String-key sorting caused the original crossed line plots. Keep numeric
        # bin edges and keep the long-read tail even when the main view is zoomed.
        counts = {"10000": 1, "500": 2, "0": 3, "1000": 0}
        points = histogram_points(counts, 500)
        self.assertEqual([x for x, n in points], [0, 500, 1000, 10000])
        self.assertEqual(sum(n for x, n in points), 6)
        r = Report(self.args)
        r.add_histogram("length", "Read-length distribution", {"batch_one": counts}, 500,
                        "Read length (nt)", "Test", focus_upper=10000)
        section = r.data["tcr_03_length"]
        self.assertEqual(section["plot_type"], "html")
        self.assertIn("5 reads below 10,000", section["data"])
        self.assertIn("1 reads at or above 10,000", section["data"])
        self.assertEqual(section["data"].count("data:image/svg+xml;base64,"), 2)
        self.assertEqual(sum(row["reads"] for row in r.tables["tcr_03_length_bins"].values()), 6)
        self.assertEqual(len(r.plot_assets), 2)
        payload = re.search(r'<script type="application/json" class="tcr-hist-data">(.*?)</script>', section["data"]).group(1)
        spec = json.loads(payload)
        self.assertEqual(spec["points"], [[0, 3], [500, 2], [1000, 0], [10000, 1]])
        self.assertEqual((spec["bin_width"], spec["xmin"], spec["xmax"], spec["focus_upper"]), (500, 0, 10500, 10000))
        self.assertIn('data-tcr-action="log"', section["data"])
        self.assertIn('data-tcr-action="svg"', section["data"])
        self.assertIn('window.Plotly', section["data"])


if __name__ == "__main__":
    unittest.main()
