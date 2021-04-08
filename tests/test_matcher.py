from primertrim.matcher import (
    CompleteMatcher, PartialMatcher,
    partial_seqs_left, partial_seqs_right,
)

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


