from primertrim.matcher import (
    PrimerMatch,
    CompleteMatcher, PartialMatcher,
)

def test_complete_match():
    m = CompleteMatcher(["TTTTTT"], 1)
    assert m.find_match("AGATTTTTT") == PrimerMatch("Complete", 3, 0, "TTTTTT")
    assert m.find_match("AATTTGTT") == PrimerMatch("Complete", 2, 1, "TTTGTT")
    assert m.find_match("AATTGGTT") == None # Two mismatches is too much


def test_partial_match():
    m = PartialMatcher(["AAAAAA"], 4)
    assert m.find_match("AAAAAGTCGT") == PrimerMatch("Partial", 0, 0, "AAAAA")
    assert m.find_match("AAAAGTCGT") == PrimerMatch("Partial", 0, 0, "AAAA")
    assert m.find_match("GTCGAAAAA") == PrimerMatch("Partial", 4, 0, "AAAAA")
    assert m.find_match("AAAGTCGGCT") == None # Length is 3, no match
