"""Artificial-read regressions for complete-local-allele reassessment.

These exercise the production evidence and decision functions. External
alignment/polishing is mocked; no experimental reads or identifiers are embedded.
"""
import argparse
from collections import Counter
from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from test_core import REF, change, read, working
from test_reassessment import evaluate
from tcr_pod5 import cmd_polish, cmd_qc, cmd_report
from tcr_pod5_core import (Settings, VARIATION_REASON, _compare_local_sequences,
                          _local_sequence, evidence, read_fasta, screen, write_fasta)
from tcr_pod5_targets import reference_digest


def fixture(insertions=False):
    """Two rare events look common only in their restricted REF/ALT subsets.

    960 reads have the dominant A deletion, 20 match the original reference,
    and 10 each have two other alleles. Both original events have 10 ALT / 20 REF,
    but all four local alleles together are 96% GAAA (or AAA for insertion sites).
    """
    ref = REF[:147] + 'CTGAAAACTGA' + REF[158:]
    sequence = ref[:150] + ref[151:]
    original, post = [], []
    for i in range(1000):
        if i < 960:
            seq, before, after = sequence, '150M1D149M', '299M'
        elif i < 980:
            seq, before, after = ref, '300M', '150M1I149M'
        elif insertions:
            base = 'T' if i < 990 else 'C'
            seq, before, after = ref[:150]+base+ref[150:], '150M1I150M', '150M1I3M1I146M'
        elif i < 990:
            seq, before, after = ref[:149]+ref[150:], '149M1D150M', '149M1D3M1I146M'
        else:
            seq, before, after = ref[:149]+ref[151:], '149M2D149M', '149M1D149M'
        original.append(read(f'r{i:04d}', seq, cigar=before, reverse=i % 2 == 1, nm=2))
        post.append(read(f'r{i:04d}', seq, cigar=after, ref='b__A5', reverse=i % 2 == 1, nm=2))
    return ref, sequence, original, post


def initial(original, ref):
    return screen(['APC']*len(original), original, {'APC': ref}, Settings())[0]


def calls(alleles):
    """Synthetic complete local observations for independent denominator tests."""
    result = {}
    for allele, n, quality in alleles:
        for _ in range(n):
            i = len(result)
            result[str(i)] = dict(status='OK', sequence=allele, minimum_quality=quality,
                                 forward_query_interval=[100, 100+len(allele)], strand='+-'[i % 2])
    return result


class LocalAlleleTests(unittest.TestCase):
    def test_competing_alleles_share_one_denominator(self):
        obs = calls([('GAAA', 88, 40), ('AAA', 6, 40), ('GAAAA', 6, 40)])
        result = _compare_local_sequences(obs, obs, 'GAAA', Settings())
        self.assertEqual(result['callable_reads'], 100)
        self.assertEqual(result['non_polished_callable_reads'], 12)
        self.assertIn('COMPETING_LOCAL_ALLELES', result['barriers'])

    def test_quality_filter_cannot_hide_a_twenty_percent_alternative(self):
        obs = calls([('GAAA', 80, 40), ('AAA', 20, 5)])
        result = _compare_local_sequences(obs, obs, 'GAAA', Settings())
        self.assertEqual(result['callable_alleles'], {'GAAA': 80})
        self.assertEqual(result['unfiltered_alternative_upper_fraction'], .2)
        self.assertIn('UNFILTERED_ALTERNATIVE_BOUND', result['barriers'])

    def test_unobservable_reads_count_against_resolution(self):
        obs = calls([('GAAA', 89, 40), ('AAA', 2, 40)])
        obs.update({f'u{i}': {'status': 'BOUNDARY_NOT_OBSERVED'} for i in range(9)})
        result = _compare_local_sequences(obs, obs, 'GAAA', Settings())
        self.assertEqual(result['callable_reads'], 91)
        self.assertEqual(result['unfiltered_alternative_upper_fraction'], .11)
        self.assertIn('UNFILTERED_ALTERNATIVE_BOUND', result['barriers'])

    def test_threshold_equality_retains_flag(self):
        obs = calls([('GAAA', 90, 40), ('AAA', 10, 40)])
        self.assertIn('COMPETING_LOCAL_ALLELES', _compare_local_sequences(obs, obs, 'GAAA', Settings())['barriers'])

    def test_losing_alternate_local_context_does_not_make_it_reference(self):
        before = calls([('GAAA', 80, 40), ('AAA', 20, 40)])
        after = deepcopy(before)
        for i in range(80, 100):
            after[str(i)].update(sequence='GAAA', forward_query_interval=[101, 105])
        result = _compare_local_sequences(before, after, 'GAAA', Settings())
        self.assertEqual(result['callable_reads'], 80)
        self.assertEqual(len(result['changed_local_observation_ids']), 20)
        self.assertIn('UNFILTERED_ALTERNATIVE_BOUND', result['barriers'])

    def test_insufficient_depth_or_one_strand_cannot_resolve(self):
        for n, one_strand in [(18, False), (100, True)]:
            obs = calls([('GAAA', n, 40)])
            if one_strand:
                for value in obs.values():
                    value['strand'] = '+'
            self.assertIn('INSUFFICIENT_SHARED_CALLABLE_SUPPORT',
                          _compare_local_sequences(obs, obs, 'GAAA', Settings())['barriers'])

    def test_majority_must_match_polished_allele(self):
        obs = calls([('AAA', 100, 40)])
        self.assertIn('LOCAL_ALLELE_DISAGREES_WITH_POLISHED_SEQUENCE',
                      _compare_local_sequences(obs, obs, 'GAAA', Settings())['barriers'])

    def test_complete_deletion_and_insertion_are_real_alleles(self):
        deletion = read('d', REF[:100]+REF[104:], cigar='100M4D196M')
        insertion = read('i', REF[:100]+'T'+REF[100:], cigar='100M1I200M')
        self.assertEqual(_local_sequence(deletion, 100, 104)['sequence'], '')
        self.assertEqual(_local_sequence(insertion, 100, 100)['sequence'], 'T')

    def test_local_call_requires_both_flanks_and_all_sequence(self):
        self.assertEqual(_local_sequence(read('r'), 0, 5)['status'], 'BOUNDARY_NOT_OBSERVED')
        self.assertEqual(_local_sequence(read('r', cigar='100M1N200M'), 99, 102)['status'], 'UNSUPPORTED_ALIGNMENT')
        self.assertEqual(_local_sequence(read('r', REF[:100]+'N'+REF[101:]), 100, 101)['status'], 'AMBIGUOUS_BASE')
        low = read('r', qual='I'*99+'!'+'I'*200)
        self.assertEqual(_local_sequence(low, 100, 104)['minimum_quality'], 0)


class IndelReassessmentTests(unittest.TestCase):
    def setUp(self):
        self.ref, self.seq, self.original, self.post = fixture()
        self.pre = initial(self.original, self.ref)

    def test_rare_adjacent_deletions_resolve_without_rewriting_original_counts(self):
        saved = deepcopy(self.pre)
        self.assertEqual(self.pre['reasons'], [VARIATION_REASON])
        sites = self.pre['variation_evidence']['sites']
        self.assertEqual({s['sequence'] for s in sites}, {'G', 'GA'})
        self.assertEqual([s['alleles'] for s in sites], [{'REF': 20, 'ALT': 10}]*2)
        self.assertEqual(len(self.pre['variation_evidence']['indel_regions']), 1)
        result = evaluate(self.pre, self.post, self.seq)
        self.assertEqual(result['status'], 'RESOLVED')
        self.assertEqual([s['original_alleles'] for s in result['sites']], [s['alleles'] for s in sites])
        counts = result['regions'][0]['local_alleles']
        self.assertEqual(counts['callable_alleles'], {'GAAA': 960, 'GAAAA': 20, 'AAAA': 10, 'AAA': 10})
        self.assertEqual(counts['non_polished_callable_fraction'], .04)
        self.assertEqual(self.pre, saved)

    def test_insertions_use_complete_competing_alleles_too(self):
        ref, sequence, original, post = fixture(insertions=True)
        pre = initial(original, ref)
        self.assertEqual([s['kind'] for s in pre['variation_evidence']['sites']], ['I', 'I'])
        result = evaluate(pre, post, sequence)
        self.assertEqual(result['status'], 'RESOLVED')
        self.assertEqual(result['regions'][0]['local_alleles']['non_polished_callable_fraction'], .04)

    def test_legacy_indel_evidence_and_incomplete_evidence_stay_review(self):
        for what in ['legacy', 'event', 'support', 'observations', 'flanks']:
            with self.subTest(what=what):
                pre = deepcopy(self.pre)
                original = pre['variation_evidence']
                if what == 'legacy':
                    original['schema_version'] = 1
                    original.pop('indel_regions')
                elif what == 'event':
                    original['sites'][0].pop('sequence')
                elif what == 'support':
                    original['sites'][0]['supporting_reads'].pop()
                elif what == 'flanks':
                    for s in original['sites'][0]['supporting_reads']:
                        if s['allele'] == 'ALT':
                            s.pop('forward_query_flanks')
                else:
                    original['indel_regions'][0]['observations'].pop('r0000')
                self.assertEqual(evaluate(pre, self.post, self.seq)['status'], 'RETAINED')

    def test_unobservable_original_alt_cannot_be_discarded(self):
        pre = deepcopy(self.pre)
        pre['variation_evidence']['indel_regions'][0]['observations']['r0980'] = {'status': 'BOUNDARY_NOT_OBSERVED'}
        result = evaluate(pre, self.post, self.seq)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('ORIGINAL_SUPPORT_NOT_EVALUABLE', result['sites'][0]['barriers'])

    def test_original_event_flanks_must_project_inside_local_region(self):
        pre = deepcopy(self.pre)
        for s in pre['variation_evidence']['sites'][0]['supporting_reads']:
            if s['allele'] == 'ALT':
                s['forward_query_flanks'] = [0, 1]
        result = evaluate(pre, self.post, self.seq)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('ORIGINAL_SUPPORT_NOT_EVALUABLE', result['sites'][0]['barriers'])

    def test_mask_or_weak_local_depth_retains_historical_flag(self):
        result = evaluate(self.pre, self.post, self.seq, extra_mask=[149])
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('MASKED_AUDITED_REGION', result['regions'][0]['barriers'])
        # Support threshold is unchanged across capture and reassessment.
        cfg = replace(Settings(), min_depth=1001)
        pre = screen(['APC']*1000, self.original, {'APC': self.ref}, cfg)[0]
        result = evaluate(pre, self.post, self.seq, cfg=cfg)
        self.assertIn('INSUFFICIENT_LOCAL_SUPPORT', result['regions'][0]['barriers'])

    def test_changed_filtered_or_missing_read_retains_flag(self):
        for post in [self.post[1:], [replace(self.post[0], mapq=0)]+self.post[1:],
                     [replace(self.post[0], qual='H'+self.post[0].qual[1:])]+self.post[1:]]:
            self.assertEqual(evaluate(self.pre, post, self.seq)['status'], 'RETAINED')

    def test_subthreshold_split_reads_still_cannot_clear_indel_flag(self):
        reads = [replace(r, tags=dict(r.tags, SA='b__A5,1,+,100M199S,60,0;')) if i < 3 else r
                 for i, r in enumerate(self.post)]
        result = evaluate(self.pre, reads, self.seq)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('ORIGINAL_READS_NOT_ALL_ELIGIBLE', result['barriers'])

    def test_structural_indel_capture_keeps_original_counts_but_disables_resolution(self):
        # Non-repeat 50-base deletion: structural, even if a polished sequence
        # eventually agrees with one population.
        reads = [read(str(i), REF[:100]+REF[150:] if i < 20 else REF,
                      cigar='100M50D150M' if i < 20 else '300M', nm=50 if i < 20 else 0,
                      reverse=i % 2 == 1) for i in range(100)]
        pre = initial(reads, REF)
        sites = pre['variation_evidence']['sites']
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0]['alleles'], {'ALT': 20, 'REF': 80})
        self.assertIn('STRUCTURAL_INDEL_NOT_REASSESSED', sites[0]['issues'])
        self.assertFalse(pre['variation_evidence']['indel_regions'])

    def test_schema_one_snp_remains_supported(self):
        from test_reassessment import deletion_fixture
        ref, seq, original, post = deletion_fixture()
        pre = initial(original, ref)
        pre['variation_evidence']['schema_version'] = 1
        pre['variation_evidence'].pop('indel_regions')
        self.assertEqual(evaluate(pre, post, seq)['status'], 'RESOLVED')

    def test_nonunique_anchors_or_structural_indel_do_not_resolve(self):
        for what in ['duplicate', 'structural']:
            pre = deepcopy(self.pre)
            if what == 'duplicate':
                pre['variation_evidence']['indel_regions'][0]['left_anchor'] = 'A'
            else:
                site = pre['variation_evidence']['sites'][0]
                site['issues'] = ['STRUCTURAL_INDEL_NOT_REASSESSED']
                site.pop('region_id')
            self.assertEqual(evaluate(pre, self.post, self.seq)['status'], 'RETAINED')

    def test_repeated_sequence_outside_target_does_not_break_unique_mapping(self):
        # Like two repeated vector segments, only one inside the inferred target.
        ref = self.ref+self.ref
        original = [replace(r, start=r.start+300) for r in self.original]
        pre = screen(['APC']*1000, original, {'APC': ref}, Settings(), targets={'APC': (300, 600)})[0]
        self.assertEqual(evaluate(pre, self.post, self.seq)['status'], 'RESOLVED')

    def test_persistent_snp_and_indel_controls_stay_flagged(self):
        # Add a balanced 12% SNP to the dominant reads, away from the rare indels.
        original = [replace(r, seq=change(r.seq, [70])) if 100 <= i < 220 else r for i, r in enumerate(self.original)]
        post = [replace(r, seq=change(r.seq, [70])) if 100 <= i < 220 else r for i, r in enumerate(self.post)]
        pre = initial(original, self.ref)
        result = evaluate(pre, post, self.seq)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('POSTPOLISH_BLOCKING_VARIATION', result['barriers'])
        # A real 80:20 one-base deletion population is still post-polish variation.
        original, post = [], []
        for i in range(100):
            s = self.ref if i < 20 else self.seq
            original.append(read(str(i), s, cigar='300M' if i < 20 else '150M1D149M', reverse=i % 2 == 1, nm=1))
            post.append(read(str(i), s, cigar='150M1I149M' if i < 20 else '299M', ref='b__A5', reverse=i % 2 == 1, nm=1))
        result = evaluate(initial(original, self.ref), post, self.seq)
        self.assertEqual(result['status'], 'RETAINED')
        self.assertIn('POSTPOLISH_BLOCKING_VARIATION', result['barriers'])


class IndelDriverTests(unittest.TestCase):
    def run_driver(self, other_reason=None, mask=False):
        ref, sequence, original, reads = fixture()
        write_fasta('refs.fasta', {'APC': ref})
        Path('settings.json').write_text('{}')
        Path('reads.sam').write_text(''.join(r.line for r in original))
        Path('trim.json').write_text(json.dumps(dict(input_reads=1000, passed_reads=1000, failed_reads=0)))
        Path('targets.json').write_text(json.dumps(dict(batch='b', reference_digest=reference_digest({'APC': ref}),
            expected_constructs=['APC'], references={'APC': dict(status='INFERRED', start=0, end=len(ref))})))
        cmd_qc(argparse.Namespace(settings='settings.json', references='refs.fasta', bam='reads.sam', samtools='samtools',
                                 inferred_targets='targets.json', trim_stats='trim.json', sample='b__A5', batch='b', well='A5'))
        pre = json.loads(Path('pass/decision.json').read_text())
        if other_reason:
            pre['reasons'].append(other_reason)
            Path('pass/decision.json').write_text(json.dumps(pre))
        if mask:
            sequence = sequence[:30]+'N'+sequence[31:]
        write_fasta('polished.fasta', {'b__A5': sequence})
        Path('dorado_polish.log').write_text('mocked integration test\n')
        with patch('tcr_pod5.dorado_version'), patch('tcr_pod5.run'), patch('tcr_pod5.do_align'), \
             patch('tcr_pod5.check_basecall_models', return_value='test-model'), \
             patch('tcr_pod5.alignments', side_effect=lambda *args: iter(reads)):
            cmd_polish(argparse.Namespace(settings='settings.json', dorado='dorado', version='2.1.2', passed='pass',
                       samtools='samtools', bam='reads.sam', threads=1, models='models', sample='b__A5', device='cpu'))
        return json.loads(Path('consensus_result/decision.json').read_text()), sequence

    def test_acceptance_warning_history_mask_and_report(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            post, sequence = self.run_driver(mask=True)
            self.assertEqual(post['status'], 'CONSENSUS_PASS')
            self.assertEqual(post['masked_bases'], 1)
            self.assertEqual(post['resolved_reasons'], [VARIATION_REASON])
            self.assertEqual(post['pre_screen_reasons'], [VARIATION_REASON])
            self.assertIn('WARN_HISTORICAL_INDEL_REASSESSED', post['warnings'])
            self.assertEqual(read_fasta('consensus_result/consensus.fasta'), {'b__A5': sequence})
            Path('demux').mkdir()
            Path('demux/manifest.json').write_text(json.dumps([dict(batch='b', well='A5', reads=1000)]))
            Path('manifest.json').write_text(json.dumps(dict(demux_dirs=['demux'], qc_dirs=['qc'], polish_dirs=['consensus_result'])))
            cmd_report(argparse.Namespace(manifest='manifest.json'))
            self.assertEqual(read_fasta('accepted_consensus.fasta'), {'b__A5': sequence})
            self.assertEqual(Path('diagnostic_consensus.fasta').read_text(), '')
            self.assertIn('WARN_HISTORICAL_INDEL_REASSESSED', Path('well_summary.tsv').read_text())

    def test_other_review_reason_survives_resolution(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            post, _ = self.run_driver(other_reason='REVIEW_UNEXPECTED_CONSTRUCT')
            self.assertEqual(post['status'], 'REVIEW')
            self.assertEqual(post['reasons'], ['REVIEW_UNEXPECTED_CONSTRUCT'])
            self.assertFalse(Path('consensus_result/consensus.fasta').exists())


if __name__ == '__main__':
    unittest.main()
