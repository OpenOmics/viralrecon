"""Assignment, masking and acceptance regressions using deterministic synthetic SAMs."""
import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from test_core import REF, change, population, read, working
from tcr_pod5 import alignments, cmd_polish, cmd_qc, cmd_report
from tcr_pod5_assignment import InformativePanel, evaluate_assignment
from tcr_pod5_core import (Settings, alignment_identity, assign, diagnostic_consensus,
                           evidence, read_fasta, screen, write_fasta)
from tcr_pod5_targets import reference_digest

_rng = random.Random(109)
LONG_REF = ''.join(_rng.choice('ACGT') for _ in range(3624))
DELETED = LONG_REF[:1000] + LONG_REF[1500:]


class DeletionTests(unittest.TestCase):
    def test_long_internal_deletion_preserves_assignment_and_sequence(self):
        r = read('del', DELETED, cigar='1000M500D2124M', nm=500)
        self.assertEqual(assign([r], {'APC': LONG_REF}, Settings(), deletion_aware=False)[0], 'low_identity')
        self.assertEqual(assign([r], {'APC': LONG_REF}, Settings())[0], 'APC')
        self.assertAlmostEqual(alignment_identity(r, Settings())['raw_identity'], 3124 / 3624)
        self.assertEqual(r.seq, DELETED)
        reads = [read(f'd{i}', DELETED, cigar=r.cigar, nm=500, reverse=i % 2 == 1) for i in range(60)]
        report, ev, _ = screen(['APC'] * 60, reads, {'APC': LONG_REF}, Settings())
        self.assertEqual(report['status'], 'PASS_SCREEN')
        self.assertEqual(ev['draft'], DELETED)
        self.assertEqual(ev['callable_fraction'], 1)

    def test_other_errors_still_count(self):
        for nm in (850, 1000):
            with self.subTest(nm=nm):
                r = read('del', DELETED, cigar='1000M500D2124M', nm=nm)
                self.assertEqual(assign([r], {'APC': LONG_REF}, Settings())[0], 'low_identity')
        short = read('short', REF[:100] + REF[110:], cigar='100M10D190M', nm=10)
        self.assertEqual(alignment_identity(short, Settings())['exempt_deleted_bases'], 0)
        insertion = read('ins', REF[:100] + 'A' * 50 + REF[100:], cigar='100M50I200M', nm=50)
        self.assertEqual(assign([insertion], {'APC': REF}, Settings())[0], 'low_identity')

    def test_terminal_or_poorly_anchored_deletion_not_exempt(self):
        for cigar, seq in [('500D3124M', LONG_REF[500:]),
                           ('10M500D3114M', LONG_REF[:10] + LONG_REF[510:])]:
            with self.subTest(cigar=cigar):
                r = read('del', seq, cigar=cigar, nm=500)
                self.assertEqual(assign([r], {'APC': LONG_REF}, Settings())[0], 'low_identity')

    def test_deletion_does_not_bypass_competitive_or_split_checks(self):
        r = read('del', DELETED, cigar='1000M500D2124M', nm=500)
        alternate = read('del', DELETED, ref='OTHER', flag=256, score=495, cigar=r.cigar, nm=500)
        self.assertEqual(assign([r, alternate], {'APC': LONG_REF, 'OTHER': LONG_REF}, Settings())[0], 'ambiguous_reference')
        r.tags['SA'] = 'OTHER,1,+,100M,60,0;'
        self.assertEqual(assign([r], {'APC': LONG_REF}, Settings())[0], 'split_alignment')

    def test_missing_or_inconsistent_nm_is_not_perfect_identity(self):
        r = read('del', DELETED, cigar='1000M500D2124M', nm=500)
        del r.tags['NM']
        self.assertEqual(assign([r], {'APC': LONG_REF}, Settings())[0], 'low_identity')
        r.tags['NM'] = '0'
        with self.assertRaisesRegex(ValueError, 'NM/CIGAR'):
            assign([r], {'APC': LONG_REF}, Settings())

    def test_intact_and_deleted_populations_remain_visible(self):
        reads = [read(f'd{i}', DELETED if i < 10 else LONG_REF,
                      cigar='1000M500D2124M' if i < 10 else None,
                      nm=500 if i < 10 else 0, reverse=i % 2 == 1) for i in range(60)]
        assigned = [assign([r], {'APC': LONG_REF}, Settings())[1] for r in reads]
        report, ev, _ = screen(['APC'] * 60, assigned, {'APC': LONG_REF}, Settings())
        self.assertEqual(report['status'], 'REVIEW')
        self.assertTrue(any(s['kind'] == 'D' and s['strong_minor'] for s in ev['sites']))

    def test_equivalent_repeat_deletions_have_coherent_draft_coordinates(self):
        ref = 'G' * 30 + 'CAAAAACT' + 'G' * 30
        # The same deleted A is represented at either edge of the homopolymer.
        reads = [read(f'r{i}', ref[:32] + ref[33:], cigar='32M1D35M' if i % 2 else '35M1D32M',
                      reverse=i % 2 == 1, nm=1) for i in range(60)]
        ev = evidence(reads, ref, Settings())
        self.assertEqual(ev['draft'], ref[:32] + ref[33:])
        self.assertEqual(ev['callable_fraction'], 1)
        self.assertFalse(ev['conflicting_fixed_indels'])

    def test_large_deletion_mixture_is_not_downgraded_by_a_homopolymer_at_its_boundary(self):
        ref = LONG_REF[:995] + 'A' * 10 + LONG_REF[1005:]
        deleted = ref[:1000] + ref[1500:]
        reads = [read(f'r{i}', deleted if i < 10 else ref, reverse=i % 2 == 1,
                      cigar='1000M500D2124M' if i < 10 else None, nm=500 if i < 10 else 0) for i in range(60)]
        report, ev, _ = screen(['APC'] * 60, reads, {'APC': ref}, Settings())
        self.assertEqual(report['status'], 'REVIEW')
        self.assertTrue(any(s['kind'] == 'D' and s['blocking_variation'] for s in ev['sites']))


class InformativeTests(unittest.TestCase):
    def panel(self):
        refs = {'APC': REF, 'OTHER': change(REF, range(60, 160))}
        return refs, InformativePanel(refs, refs)

    def test_shared_backbone_does_not_dilute_assignment(self):
        refs, panel = self.panel()
        r = read('r', refs['OTHER'], mapq=0)
        self.assertEqual(assign([r], refs, Settings())[0], 'ambiguous_reference')
        score = panel.score([r], Settings())
        self.assertEqual(score['assignment'], 'OTHER')
        self.assertEqual(score['informative_positions'], 100)
        self.assertEqual(score['margin'], 1)
        self.assertEqual(r.seq, refs['OTHER'])

    def test_backbone_only_and_deleted_markers_are_unresolved(self):
        _, panel = self.panel()
        backbone = read('r', REF[:50])
        self.assertEqual(panel.score([backbone], Settings())['assignment'], 'insufficient_informative_positions')
        deleted = read('r', REF[:60] + REF[160:], cigar='60M100D140M', nm=100)
        result = panel.score([deleted], Settings())
        self.assertEqual(result['assignment'], 'insufficient_informative_positions')
        self.assertEqual(result['deleted_informative_positions'], 100)

    def test_margin_and_identity_are_distinct_gates(self):
        refs, panel = self.panel()
        mixed = change(REF, range(60, 110))
        self.assertEqual(panel.score([read('r', mixed)], replace(Settings(), informative_min_identity=.4))['assignment'],
                         'ambiguous_informative_reference')
        bad = list(REF)
        for p in range(60, 160):
            bad[p] = next(b for b in 'ACGT' if b not in (refs['APC'][p], refs['OTHER'][p]))
        self.assertEqual(panel.score([read('r', ''.join(bad))], Settings())['assignment'], 'low_informative_identity')

    def test_low_quality_and_n_bases_do_not_add_informative_support(self):
        _, panel = self.panel()
        q = 'I' * 60 + '!' * 60 + 'I' * 180
        result = panel.score([read('r', qual=q)], replace(Settings(), min_read_q=0))
        self.assertEqual(result['informative_positions'], 40)
        self.assertEqual(result['assignment'], 'insufficient_informative_positions')
        seq = REF[:60] + 'N' * 60 + REF[120:]
        self.assertEqual(panel.score([read('r', seq)], Settings())['informative_positions'], 40)

    def test_gapped_columns_and_wrong_panels(self):
        refs, _ = self.panel()
        aligned = {'APC': REF[:60] + '-' + REF[60:], 'OTHER': refs['OTHER'][:60] + 'A' + refs['OTHER'][60:]}
        new_refs = {n: s.replace('-', '') for n, s in aligned.items()}
        panel = InformativePanel(aligned, new_refs)
        self.assertEqual(panel.excluded_variable_columns, [60])
        self.assertEqual(len(panel.columns), 100)
        with self.assertRaisesRegex(ValueError, 'original sequence'):
            InformativePanel(aligned, refs)
        with self.assertRaisesRegex(ValueError, 'exactly match'):
            InformativePanel({'APC': REF}, refs)

    def test_evaluation_reports_three_methods_without_changing_production(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            refs = {'APC': LONG_REF, 'OTHER': change(LONG_REF, range(100, 200))}
            write_fasta('refs.fasta', refs)
            write_fasta('alignment.fasta', refs)
            Path('settings.json').write_text('{}')
            r = read('del', DELETED, cigar='1000M500D2124M', nm=500)
            Path('reads.sam').write_text(r.line)
            Path('manifest.json').write_text(json.dumps(dict(batch='b', alignments=[dict(well='A1', bam='reads.sam')])))
            evaluate_assignment(argparse.Namespace(settings='settings.json', references='refs.fasta',
                reference_alignment='alignment.fasta', manifest='manifest.json', outdir='audit', samtools='samtools'), alignments)
            report = json.loads(Path('audit/assignment_evaluation.json').read_text())
            self.assertTrue(report['comparison_only'])
            self.assertEqual(report['wells']['A1']['deletion_changed_reads'], 1)
            with open('audit/read_assignment_comparison.tsv') as f:
                row = next(csv.DictReader(f, delimiter='\t'))
            self.assertEqual((row['legacy_assignment'], row['deletion_aware_assignment'], row['informative_assignment']),
                             ('low_identity', 'APC', 'APC'))
            self.assertFalse(Path('accepted_consensus.fasta').exists())


class DiagnosticTests(unittest.TestCase):
    def diagnostic(self, reads, labels=None, refs=None, **kwargs):
        refs = refs or {'APC': REF}
        report, _, selected = screen(labels or [r.ref for r in reads], reads, refs, Settings(), **kwargs)
        detail, seq = diagnostic_consensus(report, selected, refs, Settings())
        return report, detail, seq

    def test_low_support_draft_masks_unsupported_bases_without_relaxing_acceptance(self):
        reads = [read(f'r{i}', reverse=i % 2 == 1, qual='I' * 100 + '!' * 200) for i in range(12)]
        report, detail, seq = self.diagnostic(reads)
        self.assertEqual(report['status'], 'REVIEW')
        self.assertEqual(report['callable_fraction'], 0)
        self.assertFalse(report['polish_eligible'])
        self.assertTrue(detail['eligible'])
        self.assertEqual(seq, REF[:100] + 'N' * 200)
        self.assertEqual(report['settings']['min_depth'], 20)
        self.assertEqual(detail['min_depth'], 10)

    def test_unresolved_well_can_have_diagnostic_without_becoming_accepted(self):
        report, detail, seq = self.diagnostic(population(30), ['APC'] * 30 + ['unmapped'] * 20)
        self.assertFalse(report['polish_eligible'])
        self.assertIn('REVIEW_UNRESOLVED_READS', report['reasons'])
        self.assertEqual(seq, REF)
        self.assertFalse(detail['accepted'])

    def test_insufficient_or_one_strand_evidence_has_no_draft(self):
        for reads in (population(9), [read(f'r{i}') for i in range(30)]):
            with self.subTest(n=len(reads)):
                _, detail, seq = self.diagnostic(reads)
                self.assertFalse(detail['eligible'])
                self.assertIsNone(seq)

    def test_strong_single_site_variation_is_masked_even_at_eighty_percent_majority(self):
        report, detail, seq = self.diagnostic(population(100, alt=change(REF, [60]), n_alt=20))
        self.assertEqual(report['status'], 'REVIEW')
        self.assertTrue(detail['eligible'])
        self.assertEqual(seq, REF[:60] + 'N' + REF[61:])

    def test_rejected_mixtures_never_have_diagnostic(self):
        scenarios = [(population(60, alt=change(REF, [60, 160]), n_alt=20), {'APC': REF}),
                     (population(30) + [read(f'b{i}', ref='OTHER', reverse=i % 2 == 1) for i in range(30)],
                      {'APC': REF, 'OTHER': change(REF, [60])})]
        for reads, refs in scenarios:
            with self.subTest(refs=list(refs)):
                report, detail, seq = self.diagnostic(reads, refs=refs)
                self.assertEqual(report['status'], 'REJECT')
                self.assertFalse(detail['eligible'])
                self.assertIsNone(seq)

    def test_sentinel_and_unknown_target_have_no_draft(self):
        for kwargs in (dict(expected_constructs=['OTHER']), dict(targets={})):
            with self.subTest(kwargs=kwargs):
                _, detail, seq = self.diagnostic(population(30), **kwargs)
                self.assertFalse(detail['eligible'])
                self.assertIsNone(seq)

    def test_zero_assigned_reads_have_no_invented_dominant_construct(self):
        report, detail, seq = self.diagnostic([], ['unmapped'] * 30)
        self.assertIsNone(report['dominant_construct'])
        self.assertFalse(detail['eligible'])
        self.assertIsNone(seq)

    def test_qc_and_summary_expose_separate_diagnostic_path(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            write_fasta('refs.fasta', {'APC': REF})
            Path('settings.json').write_text('{}')
            Path('trim.json').write_text(json.dumps(dict(input_reads=12, passed_reads=12, failed_reads=0)))
            reads = population(12)
            Path('reads.sam').write_text(''.join(r.line for r in sorted(reads, key=lambda r: r.name)))
            Path('targets.json').write_text(json.dumps(dict(batch='b', reference_digest=reference_digest({'APC': REF}),
                expected_constructs=['APC'], references={'APC': dict(status='INFERRED', start=0, end=len(REF))})))
            cmd_qc(argparse.Namespace(settings='settings.json', references='refs.fasta', bam='reads.sam',
                   samtools='samtools', inferred_targets='targets.json', trim_stats='trim.json',
                   sample='b__A1', batch='b', well='A1'))
            self.assertFalse(Path('pass').exists())
            self.assertIn('accepted=false', Path('qc/diagnostic_consensus.fasta').read_text())
            Path('demux').mkdir()
            Path('demux/manifest.json').write_text(json.dumps([dict(well='A1', reads=12, batch='b')]))
            Path('report.json').write_text(json.dumps(dict(demux_dirs=['demux'], qc_dirs=['qc'], polish_dirs=[])))
            cmd_report(argparse.Namespace(manifest='report.json'))
            self.assertEqual(Path('accepted_consensus.fasta').read_text(), '')
            self.assertEqual(read_fasta('diagnostic_consensus.fasta'), {'b__A1': REF})
            self.assertEqual(Path('review_consensus.fasta').read_text(), '')
            summary = json.loads(Path('run_summary.json').read_text())
            self.assertEqual(summary['diagnostic_kinds'], {'DIAGNOSTIC_UNPOLISHED': 1})
            self.assertEqual(summary['statuses'], {'REVIEW': 1})

    def test_report_prefers_polished_diagnostic_and_suppresses_postpolish_rejection(self):
        for status in ('REVIEW', 'REJECT'):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as tmp, working(tmp):
                for d in ('demux', 'qc', 'pol'):
                    Path(d).mkdir()
                Path('demux/manifest.json').write_text(json.dumps([dict(well='A1', reads=30, batch='b')]))
                Path('qc/decision.json').write_text(json.dumps(dict(sample='b__A1', status='REVIEW',
                    reasons=['LOW_SUPPORT'], warnings=[], polish_eligible=True, dominant_construct='APC',
                    assigned_reads=30, callable_fraction=1, diagnostic=dict(eligible=True, callable_fraction=1))))
                write_fasta('qc/diagnostic_consensus.fasta', {'b__A1': change(REF, [60])})
                Path('pol/decision.json').write_text(json.dumps(dict(sample='b__A1', status=status,
                    reasons=['REJECT_POSTPOLISH_LINKED_HAPLOTYPES'] if status == 'REJECT' else ['REVIEW_POSTPOLISH_READ_SUPPORT'],
                    length=len(REF), masked_bases=0)))
                if status == 'REVIEW':
                    write_fasta('pol/review_consensus.fasta', {'b__A1': REF})
                Path('report.json').write_text(json.dumps(dict(demux_dirs=['demux'], qc_dirs=['qc'], polish_dirs=['pol'])))
                cmd_report(argparse.Namespace(manifest='report.json'))
                self.assertEqual(Path('accepted_consensus.fasta').read_text(), '')
                if status == 'REVIEW':
                    self.assertEqual(read_fasta('diagnostic_consensus.fasta'), {'b__A1': REF})
                    self.assertIn('DIAGNOSTIC_POLISHED', Path('diagnostic_consensus.fasta').read_text())
                else:
                    self.assertEqual(Path('diagnostic_consensus.fasta').read_text(), '')
                    self.assertEqual(Path('review_consensus.fasta').read_text(), '')


class PostPolishTests(unittest.TestCase):
    def polish(self, reads, pre_reasons=None):
        Path('pass').mkdir()
        Path('pass/decision.json').write_text(json.dumps(dict(polish_eligible=True, status='REVIEW', reasons=pre_reasons or [])))
        Path('settings.json').write_text('{}')
        Path('dorado_polish.log').write_text('mocked test\n')
        write_fasta('polished.fasta', {'sample': REF})
        with patch('tcr_pod5.dorado_version'), patch('tcr_pod5.run'), patch('tcr_pod5.do_align'), \
             patch('tcr_pod5.check_basecall_models', return_value='test-model'), \
             patch('tcr_pod5.alignments', return_value=iter(reads)):
            cmd_polish(argparse.Namespace(settings='settings.json', dorado='dorado', version='2.1.2',
                       passed='pass', samtools='samtools', bam='clean.bam', threads='1',
                       models='models', sample='b__A1', device='cpu'))
        return json.loads(Path('consensus_result/decision.json').read_text())

    def test_clean_accepted_consensus_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            reads = [replace(r, ref='sample') for r in population(60)]
            report = self.polish(reads)
            self.assertEqual(report['status'], 'CONSENSUS_PASS')
            self.assertEqual(read_fasta('consensus_result/consensus.fasta'), {'b__A1': REF})

    def test_loss_of_one_strand_blocks_acceptance(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            reads = [read(f'r{i}', ref='sample') for i in range(55)]
            reads += [read(f'u{i}', ref='*', flag=4) for i in range(5)]
            report = self.polish(reads)
            self.assertEqual(report['status'], 'REVIEW')
            self.assertIn('REVIEW_POSTPOLISH_STRAND_SUPPORT', report['reasons'])
            self.assertFalse(Path('consensus_result/consensus.fasta').exists())

    def test_polishing_cannot_clear_preexisting_mixture_evidence(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            reads = [replace(r, ref='sample') for r in population(60)]
            report = self.polish(reads, ['REVIEW_WITHIN_CONSTRUCT_VARIATION'])
            self.assertEqual(report['status'], 'REVIEW')
            self.assertIn('REVIEW_WITHIN_CONSTRUCT_VARIATION', report['reasons'])
            self.assertFalse(Path('consensus_result/consensus.fasta').exists())

    def test_linked_minor_population_rejects_even_when_marginal_sites_are_weak(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            reads = []
            for i in range(100):
                quality = list('I' * 300)
                if 50 <= i < 75:
                    quality[60] = '!'
                elif i >= 75:
                    quality[160] = '!'
                reads.append(read(f'r{i}', change(REF, [60, 160]) if i < 6 else REF,
                                  ref='sample', reverse=i % 2 == 1, qual=''.join(quality)))
            report = self.polish(reads)
            self.assertEqual(report['blocking_variation_sites'], 0)
            self.assertEqual(report['status'], 'REJECT')
            self.assertIn('REJECT_POSTPOLISH_LINKED_HAPLOTYPES', report['reasons'])
            self.assertFalse(Path('consensus_result/consensus.fasta').exists())
            self.assertFalse(Path('consensus_result/review_consensus.fasta').exists())


if __name__ == '__main__':
    unittest.main()
