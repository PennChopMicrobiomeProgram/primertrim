import os
from pathlib import Path

from primertrim.remove_primers import (
    main, CompleteMatcher, PartialMatcher,
    partial_seqs_left, partial_seqs_right,
)

DATA_DIR = Path(__file__).parent / "data"


def data_fp(filename):
    return str(DATA_DIR / filename)


def read_from(filepath):
    with open(filepath) as f:
        res = f.readlines()
    return res


def test_complete_match():
    m = CompleteMatcher(["TTTTTT"], 1)
    assert m.find_match("AGATTTTTT") == 3   # Exact match
    assert m.find_match("AATTTGTT") == 2    # One mismatch is OK
    assert m.find_match("AATTGGTT") == None # Two mismatches is too much


def test_partial_match():
    m = PartialMatcher(["AAAAAA"], 4)
    assert m.find_match("AAAAAGTCGT") == 0    # Length is 5, matches
    assert m.find_match("AAAAGTCGT") == 0     # Length is 4, matches
    assert m.find_match("GTCGAAAAA") == 4     # Length is 5, matches
    assert m.find_match("AAAGTCGGCT") == None # Length is 3, no match


def test_partial_seqs_left():
    assert list(partial_seqs_left("ABCDEFG", 3)) == \
        ["BCDEFG", "CDEFG", "DEFG", "EFG"]


def test_partial_seqs_right():
    assert list(partial_seqs_right("ABCDEFG", 3)) == \
        ["ABCDEF", "ABCDE", "ABCD", "ABC"]


def test_partial_no_results():
        assert list(partial_seqs_left("ABCDE", 5)) == []
        assert list(partial_seqs_right("ABCDE", 5)) == []


def test_main_script(tmp_path):
    input_fp = data_fp("Sub10003.V1.sputum.redo_R1.fastq")
    output_fp = str(tmp_path / "out.fastq")
    log_fp = str(tmp_path / "out.log")
    args = [
        "GCATCGATGAAGAACGCAGC",
        "-i", input_fp,
        "-o", output_fp,
        "--log", log_fp,
    ]
    main(args)

    expected_log_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_R1.log")
    assert read_from(log_fp) == read_from(expected_log_fp)

    expected_output_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_R1.fastq")
    assert read_from(output_fp) == read_from(expected_output_fp)


def test_main_script(tmp_path):
    input_fp = data_fp("Sub10003.V1.sputum.redo_R1.fastq")
    output_fp = str(tmp_path / "out.fastq")
    log_fp = str(tmp_path / "out.log")
    args = [
        "GCATCGATGAAGAACGCAGC",
        "-i", input_fp,
        "-o", output_fp,
        "--log", log_fp,
        "--num_mismatches", "1",
    ]
    main(args)

    expected_log_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_1mm_R1.log")
    assert read_from(log_fp) == read_from(expected_log_fp)

    expected_output_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_1mm_R1.fastq")
    assert read_from(output_fp) == read_from(expected_output_fp)
