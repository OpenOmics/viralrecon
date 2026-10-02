"""SAM header/model checks (mocked samtools, not executable integration)."""
import unittest
from unittest.mock import patch
import test_core  # Adds the source directory to sys.path.
from tcr_pod5 import check_basecall_models

class MetadataTests(unittest.TestCase):
    def test_header_inspection_does_not_create_new_program_record(self):
        h = "@RG\tID:a\tDS:basecall_model=model1 runid=r1\n@RG\tID:b\tDS:basecall_model=model1 runid=r2\n"
        with patch("tcr_pod5.capture", return_value=h) as call:
            self.assertEqual(check_basecall_models("input.bam", "samtools"), "model1")
            self.assertIn("--no-PG", call.call_args[0][0])

    def test_mixed_models_rejected(self):
        h = "@RG\tID:a\tDS:basecall_model=model1\n@RG\tID:b\tDS:basecall_model=model2\n"
        with patch("tcr_pod5.capture", return_value=h), self.assertRaises(ValueError):
            check_basecall_models("input.bam", "samtools")

    def test_missing_model_in_one_group_rejected(self):
        h = "@RG\tID:a\tDS:basecall_model=model1\n@RG\tID:b\tSM:test\n"
        with patch("tcr_pod5.capture", return_value=h), self.assertRaises(ValueError):
            check_basecall_models("input.bam", "samtools")
