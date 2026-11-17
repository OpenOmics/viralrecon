"""Optional real Dorado/samtools smoke test. No signal/model/GPU needed.

Runs synthetic index demultiplexing and trimming, not basecalling or polishing.
Skipped unless the requested external tools are installed.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from test_core import PRIMERS, F1, F2, R1, R2, REF, working
from tcr_pod5_core import build_scheme, rc, write_fasta
from tcr_pod5 import cmd_demux, cmd_trim, sam_to_bam

DORADO = os.environ.get('TCR_TEST_DORADO', 'dorado')
SAMTOOLS = os.environ.get('TCR_TEST_SAMTOOLS', 'samtools')
AVAILABLE = bool(shutil.which(DORADO) and shutil.which(SAMTOOLS))

@unittest.skipUnless(AVAILABLE, 'Dorado and samtools are not installed in this environment')
class ExternalTests(unittest.TestCase):
    def test_real_double_ended_demux_and_trim(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            Path('primers.tsv').write_text(PRIMERS)
            write_fasta('refs.fasta', {'APC': REF})
            build_scheme('primers.tsv', 'refs.fasta', 'scheme', 2)
            reads = {
                'a1_f': F1 + REF + rc(R1),
                'a1_r': R1 + rc(REF) + rc(F1),
                'a2_f': F2 + REF + rc(R2),
                'discordant': F1 + REF + rc(R2),
                'one_ended': F1 + REF,
            }
            with open('raw.sam', 'w') as f:
                f.write('@HD\tVN:1.6\tSO:unsorted\n')
                f.write('@RG\tID:synthetic\tSM:test\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0\n')
                for name, seq in reads.items():
                    f.write(f'{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{"I"*len(seq)}\tRG:Z:synthetic\n')
            sam_to_bam('raw.sam', 'raw.bam', SAMTOOLS)
            cmd_demux(argparse.Namespace(dorado=DORADO, samtools=SAMTOOLS, version='2.1.2',
                       scheme='scheme', bam='raw.bam', batch='synthetic', threads='2'))
            counts = {r['well']: r['reads'] for r in json.loads(Path('demux/manifest.json').read_text())}
            self.assertEqual(counts, dict(A1=2, A2=1, unclassified=2))
            Path('settings.json').write_text('{}')
            cmd_trim(argparse.Namespace(dorado=DORADO, samtools=SAMTOOLS, version='2.1.2',
                      primers='scheme/primers/A1.fasta', bam='demux/wells/A1.bam', well='A1',
                      kit='SQK-LSK114', settings='settings.json', threads='2'))
            self.assertEqual(json.loads(Path('trim_stats.json').read_text())['passed_reads'], 2)

if __name__ == '__main__': unittest.main()
