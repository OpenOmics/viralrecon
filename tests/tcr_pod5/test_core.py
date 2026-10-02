import argparse
import importlib.util
import json
import os
from pathlib import Path
import random
import sys
import tempfile
import unittest
from contextlib import contextmanager

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'bin'))
from tcr_pod5_core import (Sam, Settings, assign, build_scheme, edit_distance, evidence,
                          full_primer_match, normalize_event, parse_primers, rc, screen,
                          read_fasta, write_fasta)
from tcr_pod5 import barcode_well, check_trim, cmd_qc, cmd_report

F1 = 'CTACACGACTAACCTCGGGGGTCTTTCATTTGG'
F2 = 'CTACAGAGCTAACCTCGGGGGTCTTTCATTTGG'
R1 = 'TACAGCTGATCTCCACGATTCGGATGCAAACA'
R2 = 'TACACAGCATCTCCACGATTCGGATGCAAACA'
HEADER = 'Primer\tShorthand\tDirection\tSequence\tI-start\tI-end\tFrag-start\tFrag-end\tPair_with\n'
PRIMERS = HEADER + '\n'.join([
    f'p1\tfA1\tF\t{F1}\t3\t9\t2\t10\trA1',
    f'p2\tfA2\tF\t{F2}\t3\t9\t2\t10\trA2',
    f'p3\trA1\tR\t{R1}\t3\t9\t3\t3\t',
    f'p4\trA2\tR\t{R2}\t3\t9\t3\t3\t',
]) + '\n'
_rng = random.Random(349)
REF = ''.join(_rng.choice('ACGT') for _ in range(300))


def read(name, seq=REF, ref='APC', reverse=False, cigar=None, flag=0, mapq=60, start=0,
         score=500, nm=0, qual=None):
    flag |= 16 if reverse else 0
    q = qual if qual is not None else 'I' * len(seq)
    c = cigar or f'{len(seq)}M'
    line = f'{name}\t{flag}\t{ref}\t{start+1}\t{mapq}\t{c}\t*\t0\t0\t{seq}\t{q}\tAS:i:{score}\tNM:i:{nm}\n'
    return Sam.parse(line)


def change(seq, positions):
    arr = list(seq)
    for pos in positions:
        arr[pos] = next(x for x in 'ACGT' if x != arr[pos])
    return ''.join(arr)


def population(n=60, alt=None, n_alt=0, ref='APC', cigar=None):
    return [read(f'r{i}', alt if i < n_alt else REF, ref=ref, reverse=(i % 2 == 1),
                 cigar=cigar if i < n_alt else None) for i in range(n)]


@contextmanager
def working(path):
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


class PrimerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.tsv = self.root / 'primers.tsv'
        self.tsv.write_text(PRIMERS)

    def tearDown(self):
        self.temp.cleanup()

    def test_index_coordinates(self):
        pairs = parse_primers(self.tsv, 2)
        self.assertEqual(pairs[0]['forward']['index'], 'ACACGAC')
        self.assertEqual(pairs[0]['reverse']['index'], 'CAGCTGA')
        self.assertEqual(pairs[0]['forward']['prefix'], 'CT')

    def test_reject_incomplete_plate(self):
        with self.assertRaises(ValueError):
            parse_primers(self.tsv, 48)

    def test_reject_missing_mate(self):
        self.tsv.write_text('\n'.join(PRIMERS.splitlines()[:-1]) + '\n')
        with self.assertRaises(ValueError):
            parse_primers(self.tsv, 2)

    def test_reject_duplicate_barcode(self):
        self.tsv.write_text(PRIMERS.replace(F2, F1))
        with self.assertRaises(ValueError):
            parse_primers(self.tsv, 2)

    def test_reject_wrong_coordinates(self):
        self.tsv.write_text(PRIMERS.replace('\t3\t9\t', '\t0\t9\t'))
        with self.assertRaises(ValueError):
            parse_primers(self.tsv, 2)

    def test_scheme_and_full_reference_panel(self):
        import tomllib
        ref = self.root / 'ref.fasta'
        write_fasta(ref, {'a': REF, 'b': change(REF, [100])})
        out = self.root / 'scheme'
        build_scheme(self.tsv, ref, out, 2)
        config = tomllib.loads((out / 'barcodes.toml').read_text())
        self.assertEqual(config['arrangement']['mask2_front'], 'TA')
        self.assertEqual(config['arrangement']['barcode2_pattern'], 'RBC%02i')
        self.assertEqual(config['scoring']['flank_right_pad'], 0)
        self.assertEqual(read_fasta(out / 'reference.fasta')['a'], REF)
        self.assertIn('\tet:Z:primer\tsk:Z:any', (out / 'primers/A1.fasta').read_text())

    def test_reverse_complement(self):
        self.assertEqual(rc(R1), 'TGTTTGCATCCGAATCGTGGAGATCAGCTGTA')

    def test_barcode_label(self):
        self.assertEqual(barcode_well('TCR_WELLS_barcode01', {1: 'A1'}), 'A1')
        self.assertEqual(barcode_well('unclassified', {1: 'A1'}), 'unclassified')
        with self.assertRaises(ValueError):
            barcode_well('SQK_SOMETHING_barcode01', {1: 'A1'})


class TrimTests(unittest.TestCase):
    def test_query_global_match(self):
        errors, start, end = full_primer_match(F1, 'GATT' + F1 + REF[:60])
        self.assertEqual((errors, start, end), (0, 4, 4 + len(F1)))

    def test_insertion_tolerant_match(self):
        altered = F1[:12] + 'A' + F1[12:]
        errors, start, end = full_primer_match(F1, 'GG' + altered + REF[:50])
        self.assertEqual(errors, 1)
        self.assertEqual(start, 2)

    def test_forward_trim(self):
        raw = read('r', 'GG' + F1 + REF + rc(R1) + 'AA', flag=4)
        trimmed = read('r', REF, flag=4)
        self.assertEqual(check_trim(raw, trimmed, F1, R1, Settings()), ('pass', '+'))

    def test_reverse_trim(self):
        raw = read('r', 'GG' + R1 + rc(REF) + rc(F1) + 'AA', flag=4)
        trimmed = read('r', rc(REF), flag=4)
        self.assertEqual(check_trim(raw, trimmed, F1, R1, Settings()), ('pass', '-'))

    def test_residual_primer_rejected(self):
        raw = read('r', F1 + REF + rc(R1), flag=4)
        trimmed = read('r', F1[-5:] + REF, flag=4)
        self.assertNotEqual(check_trim(raw, trimmed, F1, R1, Settings())[0], 'pass')

    def test_one_ended_rejected(self):
        raw = read('r', F1 + REF, flag=4)
        self.assertNotEqual(check_trim(raw, read('r', REF, flag=4), F1, R1, Settings())[0], 'pass')

    def test_overtrim_rejected(self):
        raw = read('r', F1 + REF + rc(R1), flag=4)
        trimmed = read('r', REF[30:], flag=4)
        self.assertNotEqual(check_trim(raw, trimmed, F1, R1, Settings())[0], 'pass')


class AssignmentTests(unittest.TestCase):
    def test_confident_assignment(self):
        self.assertEqual(assign([read('r')], {'APC': REF}, Settings())[0], 'APC')

    def test_competitive_tie(self):
        a, b = read('r'), read('r', ref='PIK3CA', flag=256, score=495)
        self.assertEqual(assign([a, b], {'APC': REF, 'PIK3CA': REF}, Settings())[0], 'ambiguous_reference')

    def test_low_mapq(self):
        self.assertEqual(assign([read('r', mapq=0)], {'APC': REF}, Settings())[0], 'ambiguous_reference')

    def test_supplementary(self):
        a, b = read('r'), read('r', flag=2048)
        self.assertEqual(assign([a, b], {'APC': REF}, Settings())[0], 'split_alignment')

    def test_backbone_only_not_full_amplicon(self):
        r = read('r', REF[:100])
        self.assertEqual(assign([r], {'APC': REF}, Settings())[0], 'partial_target')

    def test_low_read_quality(self):
        r = read('r', qual='!' * len(REF))
        self.assertEqual(assign([r], {'APC': REF}, Settings())[0], 'low_read_quality')


class ScreenTests(unittest.TestCase):
    def decision(self, reads, refs=None, labels=None):
        return screen(labels or [r.ref for r in reads], reads, refs or {'APC': REF}, Settings())

    def test_clean_pass(self):
        report, ev, chosen = self.decision(population())
        self.assertEqual(report['status'], 'PASS_SCREEN')
        self.assertEqual(ev['draft'], REF)

    def test_fixed_snp_is_not_mixture(self):
        seq = change(REF, [60])
        report, ev, chosen = self.decision(population(alt=seq, n_alt=60))
        self.assertEqual(report['status'], 'PASS_SCREEN')
        self.assertEqual(ev['draft'], seq)

    def test_two_known_constructs_rejected(self):
        reads = population(30) + [read(f'b{i}', ref='PIK3CA', reverse=(i % 2 == 1)) for i in range(30)]
        report, _, _ = self.decision(reads, {'APC': REF, 'PIK3CA': REF})
        self.assertEqual(report['status'], 'REJECT')
        self.assertIn('REJECT_MULTIPLE_CONSTRUCTS', report['reasons'])

    def test_small_second_construct_withheld(self):
        reads = population() + [read(f'b{i}', ref='PIK3CA', reverse=(i % 2 == 1)) for i in range(3)]
        report, _, _ = self.decision(reads, {'APC': REF, 'PIK3CA': REF})
        self.assertEqual(report['status'], 'REVIEW')
        self.assertIn('REVIEW_SECOND_CONSTRUCT', report['reasons'])

    def test_same_construct_single_snp_withheld(self):
        report, _, _ = self.decision(population(alt=change(REF, [60]), n_alt=20))
        self.assertEqual(report['status'], 'REVIEW')
        self.assertIn('REVIEW_WITHIN_CONSTRUCT_VARIATION', report['reasons'])

    def test_same_construct_linked_haplotypes_rejected(self):
        report, _, _ = self.decision(population(alt=change(REF, [60, 160]), n_alt=20))
        self.assertEqual(report['status'], 'REJECT')
        self.assertIn('REJECT_LINKED_HAPLOTYPES', report['reasons'])

    def test_insufficient_reads(self):
        report, _, _ = self.decision(population(10))
        self.assertNotEqual(report['status'], 'PASS_SCREEN')
        self.assertIn('LOW_SUPPORT', report['reasons'])

    def test_one_strand_withheld(self):
        reads = [read(f'r{i}') for i in range(60)]
        report, _, _ = self.decision(reads)
        self.assertIn('REVIEW_STRAND_SUPPORT', report['reasons'])

    def test_partial_coverage_withheld(self):
        reads = [read(f'r{i}', REF[:150], reverse=(i % 2 == 1)) for i in range(60)]
        report, _, _ = self.decision(reads)
        self.assertIn('REVIEW_INCOMPLETE_COVERAGE', report['reasons'])

    def test_low_quality_alt_bases_do_not_define_haplotype(self):
        reads = population()
        seq = change(REF, [60, 160])
        q = list('I' * len(REF))
        q[60] = q[160] = '!'
        for i in range(20):
            reads[i] = read(f'r{i}', seq, reverse=(i % 2 == 1), qual=''.join(q))
        report, _, _ = self.decision(reads)
        self.assertEqual(report['status'], 'PASS_SCREEN')

    def test_polymorphic_deletion_withheld(self):
        alt = REF[:60] + REF[62:]
        reads = population(alt=alt, n_alt=20, cigar='60M2D238M')
        report, ev, _ = self.decision(reads)
        self.assertEqual(report['status'], 'REVIEW')
        self.assertTrue(any(s['kind'] == 'D' for s in ev['sites']))

    def test_fixed_deletion_consensus(self):
        alt = REF[:60] + REF[62:]
        report, ev, _ = self.decision(population(alt=alt, n_alt=60, cigar='60M2D238M'))
        self.assertEqual(report['status'], 'PASS_SCREEN')
        self.assertEqual(ev['draft'], alt)

    def test_fixed_insertion_consensus(self):
        alt = REF[:60] + 'GG' + REF[60:]
        report, ev, _ = self.decision(population(alt=alt, n_alt=60, cigar='60M2I240M'))
        self.assertEqual(report['status'], 'PASS_SCREEN')
        self.assertEqual(ev['draft'], alt)

    def test_mixture_of_indel_and_snp(self):
        seq = change(REF, [160])
        alt = seq[:60] + 'GG' + seq[60:]
        report, ev, _ = self.decision(population(alt=alt, n_alt=20, cigar='60M2I240M'))
        self.assertEqual(report['status'], 'REJECT')
        self.assertIn('REJECT_LINKED_HAPLOTYPES', report['reasons'])

    def test_unresolved_fraction_blocks_consensus(self):
        reads = population()
        report, _, _ = self.decision(reads, labels=['APC'] * 60 + ['ambiguous_reference'] * 20)
        self.assertIn('REVIEW_UNRESOLVED_READS', report['reasons'])

    def test_trim_failure_blocks_consensus(self):
        report, _, _ = screen(['APC'] * 60, population(), {'APC': REF}, Settings(),
                             dict(input_reads=100, failed_reads=40))
        self.assertIn('REVIEW_PRIMER_TRIMMING', report['reasons'])

    def test_settings_reject_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'settings.json'
            p.write_text('{"typo": 4}')
            with self.assertRaises(ValueError):
                Settings.from_json(p)

    def test_normalize_repeat_indel(self):
        self.assertEqual(normalize_event('D', 3, 'A', 'CAAAAACT'), ('D', 1, 'A'))


class DriverTests(unittest.TestCase):
    def test_qc_emits_gate_only_for_clean_wells(self):
        for mixed in (False, True):
            with self.subTest(mixed=mixed), tempfile.TemporaryDirectory() as tmp, working(tmp):
                write_fasta('refs.fasta', {'APC': REF})
                Path('settings.json').write_text('{}')
                Path('trim.json').write_text(json.dumps(dict(input_reads=60, passed_reads=60, failed_reads=0)))
                reads = population(alt=change(REF, [60, 160]), n_alt=20 if mixed else 0)
                Path('reads.sam').write_text('@HD\tVN:1.6\n' + ''.join(r.line for r in sorted(reads, key=lambda r: r.name)))
                from tcr_pod5_targets import reference_digest
                Path('targets.json').write_text(json.dumps(dict(batch='b', reference_digest=reference_digest({'APC': REF}),
                    expected_constructs=['APC'], references={'APC': dict(status='INFERRED', start=0, end=len(REF))})))
                cmd_qc(argparse.Namespace(settings='settings.json', references='refs.fasta',
                       bam='reads.sam', samtools='samtools', inferred_targets='targets.json', trim_stats='trim.json', sample='b__A1', batch='b', well='A1'))
                self.assertEqual(Path('pass').exists(), not mixed)
                self.assertEqual(Path('qc/decision.json').exists(), True)

    def test_all_rejected_report_has_no_consensus(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            Path('demux').mkdir()
            Path('qc').mkdir()
            Path('demux/manifest.json').write_text(json.dumps([
                dict(well='A1', reads=60, batch='b'), dict(well='A2', reads=0, batch='b')]))
            Path('qc/decision.json').write_text(json.dumps(dict(sample='b__A1', status='REJECT',
                 reasons=['REJECT_MULTIPLE_CONSTRUCTS'], dominant_construct='APC', assigned_reads=60, callable_fraction=1)))
            Path('report.json').write_text(json.dumps(dict(demux_dirs=['demux'], qc_dirs=['qc'], polish_dirs=[])))
            cmd_report(argparse.Namespace(manifest='report.json'))
            self.assertEqual(Path('accepted_consensus.fasta').read_text(), '')
            out = json.loads(Path('run_summary.json').read_text())
            self.assertEqual(out['accepted_consensuses'], 0)
            self.assertEqual(out['statuses']['NO_READS'], 1)



if __name__ == '__main__':
    unittest.main()
