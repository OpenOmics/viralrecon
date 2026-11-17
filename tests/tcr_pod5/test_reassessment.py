"""Reassessment regressions using artificial sequences and explicit alignments.

External Dorado/samtools calls are mocked in driver tests. No patient sequences
or reads are embedded here; these tests do not validate a polishing model.
"""
import argparse
from collections import defaultdict
from copy import deepcopy
import csv
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from test_core import REF, change, population, read, working
from tcr_pod5 import cmd_polish, cmd_qc, cmd_report
from tcr_pod5_core import (Settings, VARIATION_REASON, _read_cohort, assign,
                           evidence, rc, read_fasta, reassess_variation, screen, write_fasta)
from tcr_pod5_targets import reference_digest


def deletion_fixture():
    """A small nucleotide subset at a dominant deletion gives a false SNP flag.

    The realignment preserves every base/quality: six A calls become insertions
    and six G calls become matches to the neighboring G. These are synthetic
    alignment representations, not simulated Dorado model outputs.
    """
    ref = REF[:149] + 'GAAAAT' + REF[155:]
    consensus = ref[:150] + ref[151:]
    original, realigned = [], []
    for i in range(400):
        if i < 6:
            seq, before, after, nm = ref, '300M', '150M1I149M', 0
        elif i < 12:
            seq, before, after, nm = ref[:150] + 'G' + ref[151:], '300M', '149M1I150M', 1
        else:
            seq, before, after, nm = consensus, '150M1D149M', '299M', 1
        rid = f'r{i:03d}'
        original.append(read(rid, seq, cigar=before, nm=nm, reverse=i % 2 == 1))
        realigned.append(read(rid, seq, ref='b__A5', cigar=after, nm=1 if i < 12 else 0, reverse=i % 2 == 1))
    return ref, consensus, original, realigned


def evaluate(pre, reads, sequence, cfg=None, extra_mask=()):
    cfg = cfg or Settings()
    groups = defaultdict(list)
    for r in reads:
        groups[r.name].append(r)
    selected = []
    for records in groups.values():
        _, chosen = assign(records, {'b__A5': sequence}, cfg)
        if chosen:
            selected.append(chosen)
    ev = evidence(selected, sequence, cfg, mask_variation=True)
    masks = set(ev['masked_variation_positions']) | set(extra_mask)
    masks.update(p for p, b in enumerate(sequence) if b == 'N')
    return reassess_variation(pre, groups, selected, sequence, ev, masks, cfg)


class ReassessmentTests(unittest.TestCase):
    def setUp(self):
        self.ref, self.sequence, self.original, self.reads = deletion_fixture()
        self.pre, _, _ = screen(['APC'] * len(self.original), self.original, {'APC': self.ref}, Settings())
        self.assertEqual(self.pre['reasons'], [VARIATION_REASON])
        self.assertEqual(self.pre['blocking_variation_sites'], 1)

    def test_same_reads_resolve_deletion_alignment_artifact(self):
        record = self.pre['variation_evidence']
        self.assertEqual(record['cohort']['reads'], 400)
        self.assertEqual(record['sites'][0]['alleles'], {'A': 6, 'G': 6})
        result = evaluate(self.pre, self.reads, self.sequence)
        self.assertEqual(result['status'], 'RESOLVED')
        self.assertEqual(result['sites'][0]['projection_counts'], {'INSERTION': 6, 'MATCH': 6})
        self.assertFalse(result['barriers'])

    def test_missing_read_or_changed_id_does_not_clear_flag(self):
        for reads in (self.reads[1:], [replace(self.reads[0], name='new-id')] + self.reads[1:]):
            with self.subTest(ids=[r.name for r in reads[:1]]):
                result = evaluate(self.pre, reads, self.sequence)
                self.assertEqual(result['status'], 'RETAINED')
                self.assertIn('READ_COHORT_OR_BASECALLS_CHANGED', result['barriers'])

    def test_filtered_supporting_read_does_not_clear_flag(self):
        reads = [replace(self.reads[0], mapq=0)] + self.reads[1:]
        result = evaluate(self.pre, reads, self.sequence)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('ORIGINAL_READS_NOT_ALL_ELIGIBLE', result['barriers'])

    def test_changed_base_or_quality_does_not_clear_flag(self):
        for first in (replace(self.reads[0], seq=change(self.reads[0].seq, [20])),
                      replace(self.reads[0], qual='H' + self.reads[0].qual[1:])):
            with self.subTest(quality=first.qual[0]):
                result = evaluate(self.pre, [first] + self.reads[1:], self.sequence)
                self.assertEqual(result['status'], 'RETAINED')
                self.assertIn('READ_COHORT_OR_BASECALLS_CHANGED', result['barriers'])

    def test_cohort_fingerprint_preserves_reverse_orientation(self):
        r = self.reads[0]
        reverse = replace(r, flag=r.flag ^ 16, seq=rc(r.seq), qual=r.qual[::-1])
        self.assertEqual(_read_cohort([r]), _read_cohort([reverse]))

    def test_masked_reference_or_final_output_region_cannot_resolve(self):
        for sequence, masks in ((self.sequence[:150] + 'N' + self.sequence[151:], ()),
                                (self.sequence, (150,))):
            with self.subTest(mask=masks):
                result = evaluate(self.pre, self.reads, sequence, extra_mask=masks)
                self.assertEqual(result['status'], 'RETAINED')
                self.assertIn('MASKED_AUDITED_REGION', result['sites'][0]['barriers'])

    def test_low_local_depth_cannot_resolve(self):
        # Same configured depth before/after. A synthetic shallow observation
        # pile demonstrates the local coverage gate independently of read loss.
        groups = {r.name: [r] for r in self.reads}
        ev = evidence(self.reads, self.sequence, Settings(), mask_variation=True)
        ev['pile'][150] = dict(list(ev['pile'][150].items())[:5])
        result = reassess_variation(self.pre, groups, self.reads, self.sequence, ev, [], Settings())
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('INSUFFICIENT_LOCAL_SUPPORT', result['sites'][0]['barriers'])

    def test_changed_settings_and_missing_legacy_evidence_cannot_resolve(self):
        result = evaluate(self.pre, self.reads, self.sequence, replace(Settings(), min_depth=19))
        self.assertIn('QC_SETTINGS_CHANGED', result['barriers'])
        pre = deepcopy(self.pre)
        del pre['variation_evidence']
        result = evaluate(pre, self.reads, self.sequence)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('MISSING_OR_UNSUPPORTED_ORIGINAL_EVIDENCE', result['barriers'])

    def test_unmapped_or_nonunique_context_cannot_resolve(self):
        bad = self.sequence[:110] + change(self.sequence[110:150], [5]) + self.sequence[150:]
        result = evaluate(self.pre, self.reads, bad)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('REGION_NOT_UNIQUELY_MAPPED', result['sites'][0]['barriers'])

    def test_every_original_site_must_be_accounted_for(self):
        pre = deepcopy(self.pre)
        pre['blocking_variation_sites'] = 2
        result = evaluate(pre, self.reads, self.sequence)
        self.assertIn('INCOMPLETE_ORIGINAL_SITE_EVIDENCE', result['barriers'])
        second = deepcopy(pre['variation_evidence']['sites'][0])
        second['position_0based'] = 20
        second['issues'] = ['INSUFFICIENT_FLANKING_CONTEXT']
        pre['variation_evidence']['sites'].append(second)
        result = evaluate(pre, self.reads, self.sequence)
        self.assertEqual([s['status'] for s in result['sites']], ['RESOLVED', 'RETAINED'])
        self.assertEqual(result['status'], 'RETAINED')

    def test_unsupported_indel_and_incomplete_support_cannot_resolve(self):
        for change_record in ('indel', 'missing_support'):
            pre = deepcopy(self.pre)
            site = pre['variation_evidence']['sites'][0]
            if change_record == 'indel':
                site['kind'] = 'D'
            else:
                site['supporting_bases'].pop()
            result = evaluate(pre, self.reads, self.sequence)
            self.assertEqual(result['status'], 'RETAINED')

    def test_persistent_bidirectional_snp_controls_remain_review(self):
        # Two independent minor fractions resembling the persistent-control
        # cases; neither sequence nor count is copied from experimental data.
        for minor_count in (12, 36):
            with self.subTest(minor_count=minor_count):
                original = population(100, alt=change(REF, [160]), n_alt=minor_count)
                pre, _, _ = screen(['APC'] * 100, original, {'APC': REF}, Settings())
                reads = [replace(r, ref='b__A5') for r in original]
                result = evaluate(pre, reads, REF)
                self.assertEqual(result['status'], 'RETAINED')
                self.assertIn('POSTPOLISH_BLOCKING_VARIATION', result['barriers'])

    def test_linked_haplotypes_prevent_resolution(self):
        ev = evidence(self.reads, self.sequence, Settings(), mask_variation=True)
        ev['links'] = [{'site1': 60, 'site2': 160}]
        result = reassess_variation(self.pre, {r.name: [r] for r in self.reads}, self.reads,
                                    self.sequence, ev, [], Settings())
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('LINKED_HAPLOTYPES', result['barriers'])


class ReassessmentDriverTests(unittest.TestCase):
    def run_qc_and_polish(self, mutate=None, other_reason=False):
        ref, sequence, original, reads = deletion_fixture()
        write_fasta('refs.fasta', {'APC': ref})
        Path('settings.json').write_text('{}')
        Path('reads.sam').write_text(''.join(r.line for r in original))
        Path('trim.json').write_text(json.dumps(dict(input_reads=400, passed_reads=400, failed_reads=0)))
        Path('targets.json').write_text(json.dumps(dict(batch='b', reference_digest=reference_digest({'APC': ref}),
            expected_constructs=['APC'], references={'APC': dict(status='INFERRED', start=0, end=len(ref))})))
        cmd_qc(argparse.Namespace(settings='settings.json', references='refs.fasta', bam='reads.sam', samtools='samtools',
                                 inferred_targets='targets.json', trim_stats='trim.json', sample='b__A5', batch='b', well='A5'))
        pre = json.loads(Path('pass/decision.json').read_text())
        self.assertEqual(pre, json.loads(Path('qc/decision.json').read_text()))
        self.assertEqual(len(pre['variation_evidence']['sites'][0]['supporting_bases']), 12)
        if other_reason:
            pre['reasons'].append('REVIEW_UNEXPECTED_CONSTRUCT')
            Path('pass/decision.json').write_text(json.dumps(pre))
        if mutate:
            sequence, reads = mutate(sequence, reads)
        write_fasta('polished.fasta', {'b__A5': sequence})
        Path('dorado_polish.log').write_text('mocked driver test\n')
        with patch('tcr_pod5.dorado_version'), patch('tcr_pod5.run'), patch('tcr_pod5.do_align') as align, \
             patch('tcr_pod5.check_basecall_models', return_value='test-model'), \
             patch('tcr_pod5.alignments', side_effect=lambda *args: iter(reads)):
            cmd_polish(argparse.Namespace(settings='settings.json', dorado='dorado', version='2.1.2', passed='pass',
                       samtools='samtools', bam='reads.sam', threads=1, models='models', sample='b__A5', device='cpu'))
            self.assertEqual(align.call_count, 2)
        return json.loads(Path('consensus_result/decision.json').read_text()), sequence

    def test_qc_polish_report_accepts_resolved_flag_and_preserves_history(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            report, sequence = self.run_qc_and_polish()
            self.assertEqual(report['status'], 'CONSENSUS_PASS')
            self.assertEqual(report['resolved_reasons'], [VARIATION_REASON])
            self.assertEqual(report['pre_screen_reasons'], [VARIATION_REASON])
            self.assertEqual(report['reasons'], [])
            self.assertEqual(read_fasta('consensus_result/consensus.fasta'), {'b__A5': sequence})
            self.assertFalse(Path('consensus_result/review_consensus.fasta').exists())
            Path('demux').mkdir()
            Path('demux/manifest.json').write_text(json.dumps([dict(batch='b', well='A5', reads=400)]))
            Path('manifest.json').write_text(json.dumps(dict(demux_dirs=['demux'], qc_dirs=['qc'], polish_dirs=['consensus_result'])))
            cmd_report(argparse.Namespace(manifest='manifest.json'))
            self.assertEqual(read_fasta('accepted_consensus.fasta'), {'b__A5': sequence})
            self.assertEqual(Path('diagnostic_consensus.fasta').read_text(), '')
            with open('well_summary.tsv') as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(row['resolved_reasons'], VARIATION_REASON)
            self.assertEqual(row['variation_reassessment'], 'RESOLVED')
            self.assertEqual(json.loads(Path('run_summary.json').read_text())['resolved_review_reasons'], {VARIATION_REASON: 1})

    def test_driver_keeps_review_when_supporting_read_is_missing_or_filtered(self):
        for mutate in (lambda s, r: (s, r[1:]), lambda s, r: (s, [replace(r[0], mapq=0)] + r[1:])):
            with self.subTest(mutate=mutate), tempfile.TemporaryDirectory() as tmp, working(tmp):
                report, _ = self.run_qc_and_polish(mutate)
                self.assertEqual(report['status'], 'REVIEW')
                self.assertEqual(report['resolved_reasons'], [])
                self.assertIn(VARIATION_REASON, report['reasons'])
                self.assertFalse(Path('consensus_result/consensus.fasta').exists())

    def test_subthreshold_split_support_cannot_clear_original_variation(self):
        def split_support(sequence, reads):
            return sequence, [replace(r, tags=dict(r.tags, SA='b__A5,1,+,100M200S,60,0;'))
                              if i < 3 else r for i, r in enumerate(reads)]
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            report, _ = self.run_qc_and_polish(split_support)
            self.assertEqual(report['split_alignment_qc']['reads'], 3)
            self.assertFalse(report['split_alignment_qc']['review'])
            self.assertNotIn('REVIEW_POSTPOLISH_SPLIT_ALIGNMENTS', report['reasons'])
            self.assertEqual(report['status'], 'REVIEW')
            self.assertIn(VARIATION_REASON, report['reasons'])
            self.assertIn('ORIGINAL_READS_NOT_ALL_ELIGIBLE', report['variation_reassessment']['barriers'])
            self.assertEqual(report['resolved_reasons'], [])
            self.assertFalse(Path('consensus_result/consensus.fasta').exists())

    def test_driver_keeps_masked_site_under_review(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            report, _ = self.run_qc_and_polish(lambda s, r: (s[:150] + 'N' + s[151:], r))
            self.assertEqual(report['status'], 'REVIEW')
            self.assertIn(VARIATION_REASON, report['reasons'])
            self.assertIn('N', next(iter(read_fasta('consensus_result/review_consensus.fasta').values())))

    def test_mask_outside_audited_region_is_preserved_in_accepted_sequence(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            report, sequence = self.run_qc_and_polish(lambda s, r: (s[:30] + 'N' + s[31:], r))
            self.assertEqual(report['status'], 'CONSENSUS_PASS')
            self.assertEqual(report['masked_bases'], 1)
            self.assertEqual(read_fasta('consensus_result/consensus.fasta'), {'b__A5': sequence})

    def test_driver_persistent_snp_controls_remain_review(self):
        for minor_count in (12, 36):
            with self.subTest(minor_count=minor_count), tempfile.TemporaryDirectory() as tmp, working(tmp):
                original = population(100, alt=change(REF, [160]), n_alt=minor_count)
                pre, _, _ = screen(['APC'] * 100, original, {'APC': REF}, Settings())
                Path('pass').mkdir()
                Path('pass/decision.json').write_text(json.dumps(pre))
                Path('settings.json').write_text('{}')
                Path('dorado_polish.log').write_text('mocked persistent-variation control\n')
                write_fasta('polished.fasta', {'b__A5': REF})
                reads = [replace(r, ref='b__A5') for r in original]
                with patch('tcr_pod5.dorado_version'), patch('tcr_pod5.run'), patch('tcr_pod5.do_align'), \
                     patch('tcr_pod5.check_basecall_models', return_value='test-model'), \
                     patch('tcr_pod5.alignments', side_effect=lambda *args: iter(reads)):
                    cmd_polish(argparse.Namespace(settings='settings.json', dorado='dorado', version='2.1.2', passed='pass',
                               samtools='samtools', bam='clean.bam', threads=1, models='models', sample='b__A5', device='cpu'))
                post = json.loads(Path('consensus_result/decision.json').read_text())
                self.assertEqual(post['status'], 'REVIEW')
                self.assertEqual(post['blocking_variation_sites'], 1)
                self.assertIn('REVIEW_POSTPOLISH_VARIATION', post['reasons'])
                self.assertEqual(post['resolved_reasons'], [])

    def test_driver_preserves_other_review_reasons(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            report, _ = self.run_qc_and_polish(other_reason=True)
            self.assertEqual(report['status'], 'REVIEW')
            self.assertEqual(report['resolved_reasons'], [VARIATION_REASON])
            self.assertEqual(report['reasons'], ['REVIEW_UNEXPECTED_CONSTRUCT'])


if __name__ == '__main__':
    unittest.main()
