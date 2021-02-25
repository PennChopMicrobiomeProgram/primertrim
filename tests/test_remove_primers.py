import os
import shutil
import tempfile
import unittest

from primertrim.remove_primers import (
    main, CompleteMatcher, PartialMatcher,
    partial_seqs_left, partial_seqs_right,
)

def data_fp(filename):
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'data', filename
    )

def read_from(filepath):
    with open(filepath) as f:
        res = f.readlines()
    return res

class MatcherTests(unittest.TestCase):
    def test_1_mismatch(self):
        m = CompleteMatcher(["TTTTTT"], 1)
        self.assertEqual(m.find_match("AATTTGTT"), 2) # Actually one mismatch
        self.assertEqual(m.find_match("AGATTTTTT"), 3) # Exact match should be OK too
        self.assertEqual(m.find_match("AATTGGTT"), None) # Two mismatches is too much

    def test_partial_match(self):
        m = PartialMatcher(["AAAAAA"], 4)
        self.assertEqual(m.find_match("AAAAAGTCGT"), 0) # Length is 5, matches
        self.assertEqual(m.find_match("AAAAGTCGT"), 0) # Length is 4, still matches
        self.assertEqual(m.find_match("GTCGAAAAA"), 4) # Length is 5, matches

class FunctionTests(unittest.TestCase):
    def test_partial_seqs_left(self):
        self.assertEqual(
            list(partial_seqs_left("ABCDEFG", 3)),
            ["BCDEFG", "CDEFG", "DEFG", "EFG"])

    def test_partial_seqs_right(self):
        self.assertEqual(
            list(partial_seqs_right("ABCDEFG", 3)),
            ["ABCDEF", "ABCDE", "ABCD", "ABC"])

    def test_partial_no_results(self):
        self.assertEqual(list(partial_seqs_left("ABCDE", 5)), [])
        self.assertEqual(list(partial_seqs_right("ABCDE", 5)), [])

class ScriptTests(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.output_fp = os.path.join(self.test_dir, "out.fastq")
        self.log_fp = os.path.join(self.test_dir, "out.log")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_main_script(self):
        input_fp = data_fp("Sub10003.V1.sputum.redo_R1.fastq")
        args = [
            "GCATCGATGAAGAACGCAGC",
            "-i", input_fp,
            "-o", self.output_fp,
            "--log", self.log_fp,
        ]
        main(args)

        expected_log_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_R1.log")
        self.assertEqual(
            read_from(self.log_fp), read_from(expected_log_fp))

        expected_output_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_R1.fastq")
        self.assertEqual(
            read_from(self.output_fp), read_from(expected_output_fp))

    def test_main_script_one_mismatch(self):
        input_fp = data_fp("Sub10003.V1.sputum.redo_R1.fastq")
        args = [
            "GCATCGATGAAGAACGCAGC",
            "-i", input_fp,
            "-o", self.output_fp,
            "--log", self.log_fp,
            "--num_mismatches", "1",
        ]
        main(args)

        expected_log_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_1mm_R1.log")
        self.assertEqual(
            read_from(self.log_fp), read_from(expected_log_fp))

        expected_output_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_1mm_R1.fastq")
        self.assertEqual(
            read_from(self.output_fp), read_from(expected_output_fp))
