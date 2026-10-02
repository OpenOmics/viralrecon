"""Deterministic synthetic unit and SAM-driver integration tests; no signal data."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from test_core import REF, PRIMERS, read, change, working, population
from tcr_pod5_core import Settings, assign, evidence, screen, build_scheme, read_fasta, write_fasta
from tcr_pod5_targets import Endpoint, infer_targets, load_intervals, reference_digest, target_candidate
from tcr_pod5 import cmd_infer_targets, cmd_qc, cmd_report

ROOT = Path(__file__).resolve().parents[2]
CFG = replace(Settings(), target_min_reads=10, target_min_reads_per_well=3,
              target_min_each_strand=2, target_min_alignment_length=100, target_boundary_tolerance=10)
PANEL = {'APC': 'T'*20 + REF + 'G'*80, 'SENTINEL': 'T'*20 + change(REF, [50, 100, 150]) + 'G'*80}


def endpoints(well, n=10, start=20, end=320, ref='APC', offset=0):
    return [Endpoint(f'{well}_{offset+i}', well, ref, start, end, '+' if i%2 == 0 else '-') for i in range(n)]


def report_for(points, cfg=CFG, panel=PANEL):
    return infer_targets(panel, points, cfg, 'b', ['APC'])[0]


class InferenceTests(unittest.TestCase):
    def test_recover_interval_and_keep_sentinel(self):
        report = report_for(endpoints('A1') + endpoints('A2'))
        self.assertEqual((report['references']['APC']['start'], report['references']['APC']['end']), (20,320))
        self.assertEqual(report['references']['APC']['supporting_wells'], 2)
        self.assertEqual(report['references']['SENTINEL']['status'], 'INSUFFICIENT_SUPPORT')
        self.assertFalse(report['references']['SENTINEL']['consensus_eligible'])
        self.assertEqual(load_intervals(report, PANEL, 'b'), {'APC': (20,320)})

    def test_empty_batch(self):
        report = report_for([])
        self.assertEqual(load_intervals(report, PANEL, 'b'), {})
        self.assertTrue(all(x['start'] is None for x in report['references'].values()))

    def test_single_well_never_defines_default_target(self):
        report = report_for(endpoints('A1', 1000))
        self.assertEqual(report['references']['APC']['status'], 'INSUFFICIENT_SUPPORT')

    def test_not_enough_reads(self):
        report = report_for(endpoints('A1', 3) + endpoints('A2', 3))
        self.assertEqual(report['references']['APC']['status'], 'INSUFFICIENT_SUPPORT')

    def test_both_strands_required(self):
        points = [replace(x, strand='+') for x in endpoints('A1') + endpoints('A2')]
        self.assertEqual(report_for(points)['references']['APC']['status'], 'INSUFFICIENT_SUPPORT')

    def test_two_equal_boundary_populations_fail(self):
        points = endpoints('A1') + endpoints('A2') + endpoints('A3',start=70) + endpoints('A4',start=70)
        r = report_for(points)['references']['APC']
        self.assertEqual(r['status'], 'UNSTABLE_BOUNDARIES')
        self.assertLess(r['balanced_cluster_fraction'], .8)

    def test_depth_skew_cannot_override_multiple_wells(self):
        points = endpoints('A1', 1000, start=70)
        for i in range(2, 7):
            points += endpoints(f'A{i}', 20)
        r = report_for(points)['references']['APC']
        self.assertEqual(r['status'], 'INFERRED')
        self.assertEqual(r['start'], 20)
        self.assertEqual(r['supporting_reads'], 100)

    def test_outlier_does_not_expand_target(self):
        points = endpoints('A1', 20) + endpoints('A2', 20)
        points += endpoints('A1', 1, 0,400,offset=100) + endpoints('A2',1,0,400,offset=100)
        r=report_for(points)['references']['APC']
        self.assertEqual((r['start'], r['end']), (20,320))

    def test_tiny_well_not_used_as_boundary_seed(self):
        points = endpoints('A1') + endpoints('A2') + endpoints('A3',2,70,370)
        r=report_for(points)['references']['APC']
        self.assertEqual(r['eligible_wells'], 2)
        self.assertEqual(r['candidate_wells'], 3)

    def test_same_well_bimodality_fails(self):
        points=[]
        for w in ('A1','A2'):
            points += endpoints(w,10) + endpoints(w,10,70,370,offset=100)
        self.assertEqual(report_for(points)['references']['APC']['status'], 'UNSTABLE_BOUNDARIES')

    def test_duplicate_read_ids_fail(self):
        rows=endpoints('A1')
        with self.assertRaisesRegex(ValueError, 'Repeated read UUID'):
            report_for(rows+rows)

    def test_wrong_reference_fails(self):
        with self.assertRaisesRegex(ValueError, 'Unknown reference'):
            report_for(endpoints('A1', ref='wrong'))

    def test_out_of_range_alignment_fails(self):
        with self.assertRaisesRegex(ValueError,'outside reference'):
            report_for(endpoints('A1',end=401))

    def test_wrong_batch_cannot_reuse_target(self):
        with self.assertRaisesRegex(ValueError, 'batch/reference'):
            load_intervals(report_for([]), PANEL, 'different')

    def test_changed_sequence_cannot_reuse_target(self):
        with self.assertRaisesRegex(ValueError, 'batch/reference'):
            load_intervals(report_for([]), dict(PANEL, APC='A'*400), 'b')

    def test_missing_reference_in_report_fails(self):
        report=report_for([]);del report['references']['SENTINEL']
        with self.assertRaisesRegex(ValueError,'every reference'):
            load_intervals(report,PANEL,'b')

    def test_invalid_interval_in_report_fails(self):
        report=report_for(endpoints('A1')+endpoints('A2'))
        report['references']['APC']['end']=401
        with self.assertRaisesRegex(ValueError,'Invalid inferred'):
            load_intervals(report,PANEL,'b')

    def test_unsafe_expected_reference_fails(self):
        with self.assertRaises(ValueError):
            infer_targets(PANEL, [], CFG,'b',['missing'])

    def test_deterministic_under_input_reordering(self):
        rows=endpoints('A1')+endpoints('A2')
        a=report_for(rows); random.Random(21).shuffle(rows)
        self.assertEqual(report_for(rows), a)

    def test_realistic_reference_length(self):
        panel={'APC': 'A'*4219}
        rows=endpoints('A1',30,499,4080)+endpoints('A2',30,499,4080)
        report=infer_targets(panel,rows,Settings(),'b')[0]
        self.assertEqual(load_intervals(report,panel,'b'),{'APC':(499,4080)})


class AssignmentAndGateTests(unittest.TestCase):
    def test_full_panel_assignment_does_not_require_vector_coverage(self):
        r=read('r',start=20)
        self.assertEqual(assign([r],PANEL,CFG)[0], 'partial_target')
        self.assertEqual(assign([r],PANEL,CFG,require_target_coverage=False)[0], 'APC')

    def test_missing_mapq_255_is_not_confident(self):
        self.assertEqual(assign([read('r', mapq=255)],{'APC':REF},CFG)[0], 'ambiguous_reference')

    def test_sa_tag_is_split_evidence_even_without_supplementary_row(self):
        r=read('r');r.tags['SA']='APC,1,+,300M,60,0;'
        self.assertEqual(assign([r],{'APC':REF},CFG)[0], 'split_alignment')

    def test_secondary_score_gap_remains_enforced(self):
        r=read('r');s=read('r',ref='SENTINEL',flag=256,score=499)
        self.assertEqual(assign([r,s],PANEL,CFG,False)[0], 'ambiguous_reference')

    def test_short_read_can_assign_but_not_infer(self):
        r=read('r',seq=REF[:80])
        self.assertEqual(assign([r],PANEL,CFG,False)[0], 'APC')
        self.assertEqual(target_candidate(r,CFG), 'short_alignment')

    def test_heavily_clipped_read_cannot_define_boundary(self):
        r=read('r',cigar='20S280M')
        self.assertEqual(target_candidate(r,CFG), 'clipped_alignment')

    def test_hard_clipping_is_excluded(self):
        r=read('r',cigar='5H300M')
        self.assertEqual(assign([r],{'APC':REF},CFG)[0], 'partial_read')

    def test_callable_denominator_uses_target_not_vector(self):
        reads=[read(f'r{i}',start=20,reverse=bool(i%2)) for i in range(60)]
        report,ev,_=screen(['APC']*60,reads,PANEL,CFG,targets={'APC':(20,320)},expected_constructs=['APC'])
        self.assertEqual(report['status'],'PASS_SCREEN')
        self.assertEqual(ev['draft'],REF)
        self.assertEqual(ev['callable_fraction'],1)
        self.assertEqual(report['target_interval'],[20,320])

    def test_uninferred_target_review_no_draft(self):
        reads=population()
        report,ev,_=screen(['APC']*60,reads,{'APC':REF},CFG,targets={})
        self.assertIn('REVIEW_TARGET_NOT_INFERRED',report['reasons'])
        self.assertEqual(ev['draft'],'')
        self.assertIsNone(ev['callable_fraction'])

    def test_pure_sentinel_not_accepted(self):
        reads=population(ref='SENTINEL')
        report,_,_=screen(['SENTINEL']*60,reads,{'APC':change(REF,[100]),'SENTINEL':REF},CFG,
                            targets={'SENTINEL':(0,300)},expected_constructs=['APC'])
        self.assertEqual(report['status'],'REVIEW')
        self.assertIn('REVIEW_UNEXPECTED_CONSTRUCT',report['reasons'])

    def test_minor_sentinel_detected_without_sentinel_target(self):
        reads=population(n=50)+[read(f's{i}',ref='SENTINEL',reverse=bool(i%2)) for i in range(10)]
        report,_,_=screen(['APC']*50+['SENTINEL']*10,reads,{'APC':REF,'SENTINEL':change(REF,[100])},CFG,
                           targets={'APC':(0,300)},expected_constructs=['APC'])
        self.assertEqual(report['status'],'REJECT')
        self.assertEqual(report['construct_counts']['SENTINEL'],10)

    def test_partial_reads_do_not_redefine_target(self):
        reads=[read(f'r{i}',seq=REF[:150],reverse=bool(i%2)) for i in range(60)]
        report,_,_=screen(['APC']*60,reads,{'APC':REF},CFG,targets={'APC':(0,300)})
        self.assertEqual(report['status'],'REVIEW')
        self.assertEqual(report['target_partial_reads'],60)
        self.assertEqual(report['construct_counts']['APC'],60)

    def test_shifted_well_is_reviewed_against_batch_target(self):
        cfg=replace(CFG,target_boundary_tolerance=10)
        reads=[read(f'r{i}',seq=REF[:280],start=20,reverse=bool(i%2)) for i in range(60)]
        report,_,_=screen(['APC']*60,reads,{'APC':REF},cfg,targets={'APC':(0,300)})
        self.assertIn('REVIEW_TARGET_BOUNDARIES',report['reasons'])
        self.assertEqual(report['target_boundary_outlier_reads'],60)
        self.assertEqual(report['target_partial_reads'],0)

    def test_variants_outside_interval_not_evaluated(self):
        reads=population(alt=change(REF,[10,290]),n_alt=20)
        ev=evidence(reads,REF,CFG,interval=(20,280))
        self.assertEqual(ev['sites'],[])
        self.assertEqual(ev['draft'],REF[20:280])

    def test_variants_inside_interval_keep_original_coordinates(self):
        reads=population(alt=change(REF,[100,160]),n_alt=20)
        ev=evidence(reads,REF,CFG,interval=(20,280))
        self.assertEqual({x['pos'] for x in ev['sites']},{100,160})
        self.assertTrue(ev['links'])

    def test_nonfinite_or_invalid_config_rejected(self):
        for value in (float('nan'),float('inf'),-1,0):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as d:
                p=Path(d)/'s.json';p.write_text(json.dumps({'target_min_reads':value}))
                with self.assertRaises(ValueError): Settings.from_json(p)


class InferenceDriverTests(unittest.TestCase):
    def prepare_inputs(self, empty=False):
        write_fasta('refs.fasta',PANEL)
        Path('settings.json').write_text(json.dumps({'target_min_alignment_length':100}))
        inputs=[]
        if not empty:
            for w in ['A1','A2']:
                reads=[read(f'{w}_{i:03d}',start=20,reverse=bool(i%2)) for i in range(30)]
                Path(w+'.sam').write_text('@HD\tVN:1.6\n'+''.join(r.line for r in reads))
                inputs.append(dict(well=w,bam=w+'.sam'))
        Path('manifest.json').write_text(json.dumps(dict(batch='b',expected_constructs=['APC'],alignments=inputs)))

    def test_cli_no_bed_input_and_empty_batch_report(self):
        with tempfile.TemporaryDirectory() as d, working(d):
            self.prepare_inputs(empty=True)
            subprocess.run([sys.executable,str(ROOT/'bin/tcr_pod5.py'),'infer-targets','--batch','b',
                            '--manifest','manifest.json','--references','refs.fasta','--settings','settings.json',
                            '--outdir','targets'],check=True)
            self.assertEqual(Path('targets/inferred_targets.bed').read_text(),'')
            self.assertEqual(len(Path('targets/target_summary.tsv').read_text().splitlines()),3)

    def test_cli_inference_then_qc(self):
        with tempfile.TemporaryDirectory() as d, working(d):
            self.prepare_inputs()
            cmd_infer_targets(argparse.Namespace(settings='settings.json',references='refs.fasta',
                manifest='manifest.json',batch='b',samtools='not-needed-for-sam',outdir='targets'))
            self.assertEqual(Path('targets/inferred_targets.bed').read_text(),'APC\t20\t320\tAPC\n')
            report=json.loads(Path('targets/inferred_targets.json').read_text())
            self.assertEqual(report['input_reads'],60)
            Path('trim.json').write_text(json.dumps(dict(input_reads=30,passed_reads=30,failed_reads=0)))
            Path('settings.json').write_text(json.dumps({'min_reads':20}))
            cmd_qc(argparse.Namespace(settings='settings.json',references='refs.fasta',bam='A1.sam',
                samtools='not-needed-for-sam',trim_stats='trim.json',sample='b__A1',batch='b',well='A1',
                inferred_targets='targets/inferred_targets.json'))
            self.assertEqual(read_fasta('pass/draft.fasta'),{'b__A1':REF})

    def test_duplicate_uuid_across_wells_fails(self):
        with tempfile.TemporaryDirectory() as d, working(d):
            self.prepare_inputs()
            Path('A2.sam').write_text(Path('A1.sam').read_text())
            with self.assertRaisesRegex(ValueError,'Repeated read UUID'):
                cmd_infer_targets(argparse.Namespace(settings='settings.json',references='refs.fasta',
                    manifest='manifest.json',batch='b',samtools='samtools',outdir='targets'))

    def test_prepare_cli_generates_full_panel_without_bed(self):
        with tempfile.TemporaryDirectory() as d, working(d):
            self.prepare_inputs();Path('primers.tsv').write_text(PRIMERS)
            subprocess.run([sys.executable,str(ROOT/'bin/tcr_pod5.py'),'prepare','--primers','primers.tsv',
                '--references','refs.fasta','--outdir','scheme','--expected-pairs','2',
                '--expected-constructs','APC'],check=True,capture_output=True)
            self.assertEqual(read_fasta('scheme/reference.fasta'),PANEL)
            self.assertEqual(json.loads(Path('scheme/scheme.json').read_text())['expected_constructs'],['APC'])
            self.assertFalse(Path('scheme/amplicons.fasta').exists())

    def test_legacy_targets_cli_flag_rejected(self):
        result=subprocess.run([sys.executable,str(ROOT/'bin/tcr_pod5.py'),'prepare','--primers','p',
            '--references','r','--outdir','out','--targets','old.bed'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('unrecognized arguments',result.stderr)

    def test_all_no_read_batches_report_targets(self):
        with tempfile.TemporaryDirectory() as d, working(d):
            self.prepare_inputs(empty=True)
            cmd_infer_targets(argparse.Namespace(settings='settings.json',references='refs.fasta',
                manifest='manifest.json',batch='b',samtools='samtools',outdir='targets'))
            Path('demux').mkdir()
            Path('demux/manifest.json').write_text(json.dumps([dict(batch='b',well='A1',reads=0)]))
            Path('report.json').write_text(json.dumps(dict(demux_dirs=['demux'],qc_dirs=[],polish_dirs=[],target_dirs=['targets'])))
            cmd_report(argparse.Namespace(manifest='report.json'))
            r=json.loads(Path('run_summary.json').read_text())
            self.assertIn('b',r['batch_targets'])
            self.assertEqual(r['accepted_consensuses'],0)


if __name__ == '__main__':
    unittest.main()
