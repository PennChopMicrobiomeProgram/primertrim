from primertrim.trimmable_reads import TrimmableReads


read1 = ("seq1", "ATGTCATGACTTGACTGCGG", "FFFFFFFFFFFFFFFFFFFF")
read2 = ("seq2", "AGTCACGCTGACTGCATTGA", "FFFFFFFFFFFFFFFFFFFF")
read3 = ("seq3", "TACGTCATGCATCGTAGTAA", "FFFFFFFFFFFFFFFFFFFF")

seq1 = ("seq1", "ATGTCATGACTTGACTGCGG")
seq2 = ("seq2", "AGTCACGCTGACTGCATTGA")
seq3 = ("seq3", "TACGTCATGCATCGTAGTAA")

log1 = ("seq1", "No match", 20, None, None)
log2 = ("seq2", "No match", 20, None, None)
log3 = ("seq3", "No match", 20, None, None)

read2_trim10 = ("seq2", "AGTCACGCTG", "FFFFFFFFFF")
log2_trim10 = ("seq2", "Complete", 10, 0, "ACTGCATTGA")
read2_trim0 = ("seq2", "", "")


class MockMatch:
    method = "Complete"
    mismatches = 0
    primerseq = "ACTGCATTGA"

    def __init__(self, start):
        self.start = start


def test_register_match():
    t = TrimmableReads([read1, read2, read3])

    assert list(t.get_unmatched_seqs()) == [seq1, seq2, seq3]

    m = MockMatch(10)
    t.register_match("seq2", m)

    assert list(t.get_unmatched_seqs()) == [seq1, seq3]


def test_output():
    t = TrimmableReads([read1, read2, read3])

    assert list(t.output_reads()) == [read1, read2, read3]

    m = MockMatch(10)
    t.register_match("seq2", m)

    assert list(t.output_reads()) == [read1, read2_trim10, read3]
    assert list(t.output_loginfo()) == [log1, log2_trim10, log3]


def test_output_min_length():
    t = TrimmableReads([read1, read2, read3])
    m = MockMatch(10)
    t.register_match("seq2", m)

    # All reads written out
    assert list(t.output_reads(min_length=0)) == [read1, read2_trim10, read3]
    # Removes read2
    assert list(t.output_reads(min_length=15)) == [read1, read3]
    # No reads are long enough
    assert list(t.output_reads(min_length=30)) == []


def test_output_zero_length():
    t = TrimmableReads([read1, read2, read3])
    m = MockMatch(0)
    t.register_match("seq2", m)

    # All reads written out
    assert list(t.output_reads(min_length=0)) == [read1, read2_trim0, read3]
    # Negative min_length writes out all reads
    assert list(t.output_reads(min_length=-5)) == [read1, read2_trim0, read3]
    # Removes read2
    assert list(t.output_reads(min_length=1)) == [read1, read3]
