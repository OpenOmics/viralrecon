"""Split-read review gates and consensus exclusions through the stage drivers."""
import argparse
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from test_core import REF, population, read, working
from tcr_pod5 import cmd_polish, cmd_qc
from tcr_pod5_core import read_fasta, write_fasta
from tcr_pod5_targets import reference_digest


def split_records(n, ref='APC', kind='SA'):
    records = []
    for i in range(n):
        if kind == 'N':
            r = read(f'split{i}', REF[:100] + REF[110:], ref=ref,
                     cigar='100M10N190M', reverse=i % 2 == 1)
        else:
            r = read(f'split{i}', ref=ref, reverse=i % 2 == 1)
        if kind == 'SA':
            r.line = r.line.rstrip('\n') + f'\tSA:Z:{ref},1,+,100M200S,60,0;\n'
        records.append(r)
        if kind == 'supplementary':
            records.append(read(r.name, ref=ref, flag=2048, reverse=i % 2 == 1))
    return records


def initial_qc(records, settings=None):
    write_fasta('refs.fasta', {'APC': REF})
    Path('settings.json').write_text(json.dumps(settings or {}))
    total = len({r.name for r in records})
    Path('trim.json').write_text(json.dumps(dict(input_reads=total, passed_reads=total, failed_reads=0)))
    Path('reads.sam').write_text(''.join(r.line for r in sorted(records, key=lambda r: r.name)))
    Path('targets.json').write_text(json.dumps(dict(batch='b', reference_digest=reference_digest({'APC': REF}),
        expected_constructs=['APC'], references={'APC': dict(status='INFERRED', start=0, end=len(REF))})))
    cmd_qc(argparse.Namespace(settings='settings.json', references='refs.fasta', bam='reads.sam',
        samtools='samtools', inferred_targets='targets.json', trim_stats='trim.json',
        sample='b__A1', batch='b', well='A1'))
    return json.loads(Path('qc/decision.json').read_text())


def final_qc(records, settings=None):
    Path('pass').mkdir()
    Path('pass/decision.json').write_text(json.dumps(dict(polish_eligible=True, status='PASS_SCREEN', reasons=[])))
    Path('settings.json').write_text(json.dumps(settings or {}))
    Path('dorado_polish.log').write_text('mocked test\n')
    write_fasta('polished.fasta', {'APC': REF})
    # Round-trip the SAM text so tagged records are parsed exactly as in cmd_qc.
    from tcr_pod5_core import Sam
    with patch('tcr_pod5.dorado_version'), patch('tcr_pod5.run'), patch('tcr_pod5.do_align'), \
            patch('tcr_pod5.check_basecall_models', return_value='test-model'), \
            patch('tcr_pod5.alignments', return_value=iter([Sam.parse(r.line) for r in records])):
        cmd_polish(argparse.Namespace(settings='settings.json', dorado='dorado', version='2.1.2',
            passed='pass', samtools='samtools', bam='clean.bam', threads='1',
            models='models', sample='b__A1', device='cpu'))
    return json.loads(Path('consensus_result/decision.json').read_text())


class SplitReadTests(unittest.TestCase):
    def test_count_and_fraction_are_both_required_before_and_after_polish(self):
        # Count only, fraction only, exact boundary, below boundary, and neither.
        cases = [(997, 3, False), (58, 2, False), (97, 3, True),
                 (98, 3, False), (99, 1, False), (60, 0, False)]
        for driver, flag, passing in ((initial_qc, 'REVIEW_SPLIT_ALIGNMENTS', 'PASS_SCREEN'),
                                      (final_qc, 'REVIEW_POSTPOLISH_SPLIT_ALIGNMENTS', 'CONSENSUS_PASS')):
            for n_clean, n_split, blocked in cases:
                with self.subTest(stage=driver.__name__, clean=n_clean, split=n_split), \
                        tempfile.TemporaryDirectory() as tmp, working(tmp):
                    report = driver(population(n_clean) + split_records(n_split))
                    self.assertEqual(flag in report['reasons'], blocked)
                    self.assertEqual(report['status'], 'REVIEW' if blocked else passing)
                    self.assertEqual(report['selected_reads'], n_clean)
                    metric = report['split_alignment_qc']
                    self.assertEqual(metric['reads'], n_split)
                    self.assertEqual(metric['total_reads'], n_clean + n_split)
                    self.assertAlmostEqual(metric['fraction'], n_split / (n_clean + n_split))
                    self.assertEqual(metric['review'], blocked)
                    self.assertEqual((metric['min_reads'], metric['min_fraction']), (3, .03))
                    if driver is initial_qc:
                        self.assertEqual(report['polish_eligible'], not blocked)
                        self.assertEqual(Path('pass/accepted.ids').exists(), not blocked)
                    else:
                        self.assertEqual(Path('consensus_result/consensus.fasta').exists(), not blocked)

    def test_all_evaluated_read_ids_are_the_denominator(self):
        # Three splits / 103 total is below 3%; using only 97 panel-assigned
        # reads, or excluding the ambiguous reads, would incorrectly flag it.
        for driver in (initial_qc, final_qc):
            with self.subTest(stage=driver.__name__), tempfile.TemporaryDirectory() as tmp, working(tmp):
                records = population(97) + split_records(3)
                records += [read(f'ambiguous{i}', mapq=0) for i in range(3)]
                report = driver(records)
                self.assertEqual(report['split_alignment_qc']['total_reads'], 103)
                self.assertAlmostEqual(report['split_alignment_qc']['fraction'], 3 / 103)
                self.assertFalse(report['split_alignment_qc']['review'])
                self.assertEqual(report['ambiguous_reference_reads'], 3)
                self.assertEqual(report['unresolved_for_qc'], 3)

    def test_supplementary_records_do_not_dilute_fraction_or_inflate_read_count(self):
        for driver in (initial_qc, final_qc):
            with self.subTest(stage=driver.__name__), tempfile.TemporaryDirectory() as tmp, working(tmp):
                report = driver(population(97) + split_records(3, kind='supplementary'))
                self.assertEqual(report['split_alignment_qc']['reads'], 3)
                self.assertEqual(report['split_alignment_qc']['total_reads'], 100)
                self.assertEqual(report['split_alignment_qc']['fraction'], .03)
                self.assertTrue(report['split_alignment_qc']['review'])
                self.assertEqual(report['selected_reads'], 97)

    def test_subthreshold_splits_cannot_supply_consensus_depth(self):
        for driver in (initial_qc, final_qc):
            for kind in ('SA', 'supplementary', 'N'):
                with self.subTest(stage=driver.__name__, kind=kind), \
                        tempfile.TemporaryDirectory() as tmp, working(tmp):
                    # Nineteen valid observations at position 60: including the
                    # three excluded split reads would incorrectly make it callable.
                    clean = [read(f'r{i}', reverse=i % 2 == 1,
                                  qual='I' * 60 + ('I' if i < 19 else '!') + 'I' * 239)
                             for i in range(100)]
                    report = driver(clean + split_records(3, kind=kind))
                    self.assertFalse(report['split_alignment_qc']['review'])
                    self.assertEqual(report['selected_reads'], 100)
                    self.assertEqual(report['strand_counts'], {'+': 50, '-': 50})
                    if driver is initial_qc:
                        ids = Path('pass/accepted.ids').read_text().splitlines()
                        self.assertEqual(set(ids), {r.name for r in clean})
                        fasta = 'pass/draft.fasta'
                    else:
                        self.assertEqual(report['masked_bases'], 1)
                        fasta = 'consensus_result/consensus.fasta'
                    self.assertEqual(read_fasta(fasta), {'b__A1': REF[:60] + 'N' + REF[61:]})

    def test_configured_thresholds_are_used_at_both_stages(self):
        for driver in (initial_qc, final_qc):
            for settings, blocked in (({'alert_reads': 5}, False),
                                      ({'alert_fraction': .04}, False),
                                      ({'alert_reads': 5, 'alert_fraction': .05}, True)):
                n_split = 5 if settings.get('alert_fraction') == .05 else 3
                with self.subTest(stage=driver.__name__, settings=settings), \
                        tempfile.TemporaryDirectory() as tmp, working(tmp):
                    report = driver(population(100 - n_split) + split_records(n_split), settings)
                    self.assertEqual(report['split_alignment_qc']['review'], blocked)
                    self.assertEqual(report['split_alignment_qc']['min_reads'], settings.get('alert_reads', 3))
                    self.assertEqual(report['split_alignment_qc']['min_fraction'], settings.get('alert_fraction', .03))

    def test_unresolved_fraction_remains_a_separate_blocker(self):
        for driver, reason in ((initial_qc, 'REVIEW_UNRESOLVED_READS'),
                               (final_qc, 'REVIEW_POSTPOLISH_ALIGNMENTS')):
            with self.subTest(stage=driver.__name__), tempfile.TemporaryDirectory() as tmp, working(tmp):
                # 2% split + 9% unmapped exceeds the 10% unresolved limit.
                records = population(89) + split_records(2)
                records += [read(f'unmapped{i}', flag=4, ref='*') for i in range(9)]
                report = driver(records)
                self.assertFalse(report['split_alignment_qc']['review'])
                self.assertEqual(report['unresolved_for_qc'], 11)
                self.assertIn(reason, report['reasons'])
                self.assertEqual(report['status'], 'REVIEW')

    def test_empty_initial_input_has_zero_split_fraction(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            report = initial_qc([])
            self.assertEqual(report['split_alignment_qc']['total_reads'], 0)
            self.assertEqual(report['split_alignment_qc']['fraction'], 0)
            self.assertFalse(report['split_alignment_qc']['review'])
            self.assertFalse(report['polish_eligible'])


if __name__ == '__main__':
    unittest.main()
