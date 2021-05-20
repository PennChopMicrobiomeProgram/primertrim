from primertrim.fastq import TrimmableReads

from collections import namedtuple


input_reads = [
    ("seq1", "ATGTCATGACTTGACTGCGG", "FFFFFFFFFFFFFFFFFFFF"),
    ("seq2", "AGTCACGCTGACTGCATTGA", "FFFFFFFFFFFFFFFFFFFF"),
    ("seq3", "TACGTCATGCATCGTAGTAA", "FFFFFFFFFFFFFFFFFFFF"),
]

output_reads = [
    ("seq1", "ATGTCATGACTTGACTGCGG", "FFFFFFFFFFFFFFFFFFFF"),
    ("seq2", "AGTCACGCTG", "FFFFFFFFFF"),
    ("seq3", "TACGTCATGCATCGTAGTAA", "FFFFFFFFFFFFFFFFFFFF"),
]

output_loginfo = [
    ("seq1", "No match", 20, "", ""),
    ("seq2", "Complete", 10, 0, "ACTGCATTGA"),
    ("seq3", "No match", 20, "", ""),
]

class MockMatch:
    method = "Complete"
    start = 10
    mismatches = 0
    primerseq = "ACTGCATTGA"


def test_register_match():
    t = TrimmableReads(input_reads)

    read_ids = [read_id for read_id, seq in t.get_unmatched_seqs()]
    assert read_ids == ["seq1", "seq2", "seq3"]

    m = MockMatch()
    t.register_match("seq2", m)

    read_ids = [read_id for read_id, seq in t.get_unmatched_seqs()]
    assert read_ids == ["seq1", "seq3"]


def test_output():
    t = TrimmableReads(input_reads)

    # No matches found, nothing changed from input
    assert list(t.output_reads()) == input_reads

    m = MockMatch()
    t.register_match("seq2", m)

    # Now we have the expected output
    assert list(t.output_reads()) == output_reads
    assert list(t.output_loginfo()) == output_loginfo

