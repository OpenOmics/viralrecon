"""POD5 audit logic tests with a mocked POD5 reader; not real signal decoding."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from test_core import working
from tcr_pod5 import cmd_audit

class AuditTests(unittest.TestCase):
    def perform(self, duplicates=False, empty=False, allow=False, kit='SQK-LSK114'):
        info = SimpleNamespace(acquisition_id='run', flow_cell_id='flow', sample_id='sample',
                               sequencing_kit='SQK-LSK114', sample_rate=5000)
        class Reader:
            def __init__(self, path): self.path = path
            def __enter__(self): return self
            def __exit__(self, *args): pass
            def reads(self):
                rid = 'same-read' if duplicates else str(self.path)
                yield SimpleNamespace(read_id=rid, run_info=info)
        fake = SimpleNamespace(Reader=Reader, __version__='mocked-SDK')
        entries = []
        for i in range(2):
            path = Path(f'part{i}.pod5')
            path.write_bytes(b'' if empty and i else b'synthetic-container-placeholder')
            entries.append(dict(staged=str(path), source=str(path), batch=f'b{i}', kit=kit))
        # Empty files can be skipped only if their batch still has real reads.
        if empty:
            entries[1]['batch'] = 'b0'
        Path('input.json').write_text(json.dumps(entries))
        with patch.dict(sys.modules, {'pod5': fake}):
            cmd_audit(argparse.Namespace(manifest='input.json', allow_empty=allow))

    def test_disjoint_ids_pass(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            self.perform()
            self.assertTrue(Path('audit.ok').is_file())

    def test_duplicate_ids_fail(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            with self.assertRaises(ValueError): self.perform(duplicates=True)
            self.assertFalse(Path('audit.ok').exists())
            self.assertTrue(json.loads(Path('pod5_audit.json').read_text())['duplicate_reads'])

    def test_empty_files_fail_by_default(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            with self.assertRaises(ValueError): self.perform(empty=True)

    def test_empty_files_explicitly_ignored(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            self.perform(empty=True, allow=True)
            self.assertTrue(Path('audit.ok').exists())

    def test_kit_mismatch_fail(self):
        with tempfile.TemporaryDirectory() as tmp, working(tmp):
            with self.assertRaises(ValueError): self.perform(kit='SQK-OTHER')

if __name__ == '__main__': unittest.main()
