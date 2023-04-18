from primertrim.matcher import (
    PrimerMatch,
    CompleteMatcher,
    PartialMatcher,
    AlignmentMatcher,
)

import os


def test_complete_match():
    m = CompleteMatcher(["TTTTTT"], 1, False)
    assert m.find_match("AGATTTTTT") == PrimerMatch("Complete", 3, 0, "TTTTTT")
    assert m.find_match("AATTTGTT") == PrimerMatch("Complete", 2, 1, "TTTGTT")
    assert m.find_match("AATTGGTT") == None  # Two mismatches is too much


def test_partial_match():
    m = PartialMatcher(["AAAAAA"], 4, False)
    assert m.find_match("AAAAAGTCGT") == PrimerMatch("Partial", 0, 0, "AAAAA")
    assert m.find_match("AAAAGTCGT") == PrimerMatch("Partial", 0, 0, "AAAA")
    assert m.find_match("GTCGAAAAA") == PrimerMatch("Partial", 4, 0, "AAAAA")
    assert m.find_match("AAAGTCGGCT") == None  # Length is 3, no match


def test_align_match(tmp_path):
    align_fp = str(tmp_path / "vsearch_align_fp")
    os.mkdir(align_fp)
    test_dict = {
        "test": "GGGGGAAAAA",
        "test2": "GGGGGAACAA",
        "test3": "GGGGGGGAAAAA",
        "test4": "GGGGGGGCAAAACCCCCTTTTT",
    }
    p_match = [
        PrimerMatch("Alignment", 0, 0, "GGGGGAAAAA"),
        PrimerMatch("Alignment", 0, 1, "GGGGGAACAA"),
        PrimerMatch("Alignment", 2, 0, "GGGGGAAAAA"),
        PrimerMatch("Alignment", 2, 1, "GGGGGCAAAACCCCCTTTTT"),
    ]
    m = AlignmentMatcher(["GGGGGAAAAACCCCCTTTTT"], align_fp, 0.8)
    m_gen = m.find_in_seqs(test_dict.items())
    i = 0
    for item in m_gen:
        assert item[1] == p_match[i]
        i += 1
