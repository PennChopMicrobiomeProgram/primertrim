from primertrim.dna import (
    partial_seqs_left, partial_seqs_right,
)

def test_partial_seqs_left():
    assert list(partial_seqs_left("ABCDEFG", 3)) == \
        ["BCDEFG", "CDEFG", "DEFG", "EFG", "GG"]


def test_partial_seqs_right():
    assert list(partial_seqs_right("ABCDEFG", 3)) == \
        ["ABCDEF", "ABCDE", "ABCD", "ABC"]


def test_partial_no_results():
        assert list(partial_seqs_left("ABCDE", 5)) == []
        assert list(partial_seqs_right("ABCDE", 5)) == []
