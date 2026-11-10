"""Check staging with repeated filenames, multiple batches, and empty outputs."""
import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'assets/tcr_pod5/multiqc'))
from collect_tcr_multiqc import Report
from stage_tcr_multiqc import stage


class StagingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.manifest = self.root / 'manifest.json'
        self.out = self.root / 'collected'

    def write_manifest(self, entries):
        self.manifest.write_text(json.dumps(entries))

    def test_repeated_basenames_keep_batch_and_well_identity(self):
        entries = []
        for batch, well, n in [('one', 'A1', 10), ('one', 'A2', 20), ('two', 'A1', 30)]:
            source = self.root / batch / well / 'decision.json'
            source.parent.mkdir(parents=True)
            source.write_text(json.dumps({'selected_reads': n}))
            entries.append({'source': str(source),
                            'destination': f'{batch}/wells/{well}/heterogeneity/qc/decision.json'})
        self.write_manifest(entries)
        self.assertEqual(stage(self.manifest, self.out), 3)
        for entry, n in zip(entries, [10, 20, 30]):
            self.assertEqual(json.loads((self.out / entry['destination']).read_text())['selected_reads'], n)

    def test_missing_or_duplicate_input_fails_before_output_exists(self):
        source = self.root / 'exists.json'
        source.write_text('{}')
        entry = {'source': str(source), 'destination': 'one/decision.json'}
        for entries in ([entry, entry], [entry, {'source': str(self.root / 'missing.json'), 'destination': 'two/decision.json'}]):
            self.write_manifest(entries)
            with self.assertRaises(ValueError):
                stage(self.manifest, self.out)
            self.assertFalse(self.out.exists())

    def test_invalid_destinations_are_rejected(self):
        source = self.root / 'exists.json'
        source.write_text('{}')
        for name in ('../escape.json', '/absolute.json', 'one/../../escape', 'one//decision.json'):
            self.write_manifest([{'source': str(source), 'destination': name}])
            with self.assertRaises(ValueError):
                stage(self.manifest, self.out)
            self.assertFalse(self.out.exists())

    def test_zero_read_run_needs_no_well_or_polish_inputs(self):
        sources = self.root / 'sources'
        sources.mkdir()
        payloads = {
            'one/demultiplexing/demux/manifest.json': json.dumps([
                {'batch': 'one', 'well': 'A1', 'reads': 0},
                {'batch': 'one', 'well': 'unclassified', 'reads': 10}]),
            'one/basecalling/basecall_provenance.json': json.dumps({'batch': 'one', 'reads': 10, 'dorado_version': 'TEST'}),
            'one/target_inference/targets/inferred_targets.json': json.dumps({'batch': 'one', 'settings': {'min_reads': 50}, 'references': {}}),
            'summary/well_summary.tsv': 'sample\tbatch\twell\tstatus\traw_assigned_reads\tdiagnostic_available\n'
                                        'one__A1\tone\tA1\tNO_READS\t0\tfalse\n',
            'summary/run_summary.json': json.dumps({'wells': 1, 'statuses': {'NO_READS': 1},
                'accepted_consensuses': 0, 'review_consensuses': 0, 'diagnostic_consensuses': 0}),
            'summary/accepted_consensus.fasta': '',
            'summary/review_consensus.fasta': '',
            'summary/diagnostic_consensus.fasta': '',
            # Even if an audit exists, production report arguments must exclude it.
            'one/assignment_evaluation.json': json.dumps({'batch': 'one', 'reads': 10, 'wells': {}}),
        }
        entries = []
        for i, (destination, content) in enumerate(payloads.items()):
            source = sources / str(i)
            source.write_text(content)
            entries.append({'source': str(source), 'destination': destination})
        self.write_manifest(entries)
        stage(self.manifest, self.out)
        args = argparse.Namespace(results_dir=str(self.out), batch=None, audit_dir=[],
            assignment_evaluation=[], no_discover_audits=True)
        report = Report(args)
        report.load(); report.overview(); report.inputs_and_basecalls()
        report.demux_trim_assignment(); report.target_and_qc(); report.consensus()
        report.audits(); report.completeness()
        self.assertTrue(all(int(key.split('_')[1]) <= 11 for key in report.data))
        self.assertIn('tcr_10_files', report.data)
        self.assertIn('tcr_11_settings', report.data)
        self.assertEqual(report.data['tcr_10_files']['data']['accepted_consensus.fasta']['sequences'], 0)
        self.assertEqual(report.availability['one / A1']['qc'], 'NOT EXPECTED')
        self.assertEqual(report.data['tcr_00_batch']['data']['one']['accepted_consensuses'], 0)
        self.assertEqual(report.warnings, [])


if __name__ == '__main__':
    unittest.main()
