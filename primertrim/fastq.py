class TrimmableReads:
    def __init__(self, reads):
        self.descs = {}
        self.seqs = {}
        self.quals = {}
        self.matches = {}
        for desc, seq, qual in reads:
            read_id = get_read_id(desc)
            if read_id in self.descs:
                raise ValueError("Duplicate read ID: {}",format(read_id))
            self.descs[read_id] = desc
            self.seqs[read_id] = seq
            self.quals[read_id] = qual
            self.matches[read_id] = None

    @classmethod
    def from_fastq(cls, f):
        return cls(parse_fastq(f))

    def get_unmatched_seqs(self):
        for read_id, match in self.matches.items():
            if match is None:
                yield read_id, self.seqs[read_id]

    def register_match(self, read_id, matchobj):
        self.matches[read_id] = matchobj

    def get_trimmed_reads(self):
        for read_id in self.descs.keys():
            matchobj = self.matches[read_id]
            desc = self.descs[read_id]
            seq = self.seqs[read_id]
            qual = self.quals[read_id]
            if matchobj is not None:
                seq = seq[:matchobj.start]
                qual = qual[:matchobj.start]
            yield (desc, seq, qual)


def parse_fastq(f):
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield (desc, seq, qual)


def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return zip(*args)


def get_read_id(desc):
    return desc.split(maxsplit=1)[0]
