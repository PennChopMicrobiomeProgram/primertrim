import itertools
import sys
import argparse
import os

class FastqRead(object):
   def __init__(self, desc, seq, qual):
      self.desc = desc
      self.seq = seq
      self.qual = qual

   def trim(self, idx):
      if idx is None:
         return self
      else:
         return self.__class__(self.desc, self.seq[:idx], self.qual[:idx])

   def format_fastq(self):
      return "@{0}\n{1}\n+\n{2}\n".format(self.desc, self.seq, self.qual)

   @classmethod
   def parse(cls, f):
      for desc, seq, _, qual in _grouper(f, 4):
         desc = desc.rstrip()[1:]
         seq = seq.rstrip()
         qual = qual.rstrip()
         yield cls(desc, seq, qual)

def _grouper(iterable, n):
   "Collect data into fixed-length chunks or blocks"
   # grouper('ABCDEFG', 3) --> ABC DEF
   args = [iter(iterable)] * n
   return zip(*args)

class Matcher(object):
    def __init__(self, queryset):
        self.queryset = queryset.copy()

    def find_match(self, seq):
        idx = -1
        for query in self.queryset:
            idx = seq.find(query)
            if idx != -1:
                break
        if idx == -1:
            idx = None
        return idx


class CompleteMatcher(Matcher):
    def __init__(self, queryset, max_mismatch):
        super().__init__(queryset)

        if max_mismatch > 0:
           n = 1
           while n <= max_mismatch:
              self.queryset.extend(self._mismatched_queries(queryset, n))
              n += 1

    def _mismatched_queries(self, queryset, n_mismatch):
        # This algorithm is terrible unless the number of mismatches is very small
        assert(n_mismatch in [1, 2, 3])
        for query in queryset:
            idx_sets = itertools.combinations(range(len(query)), n_mismatch)
            for idx_set in idx_sets:
                # Replace base at each position with a literal "N", to match
                # ambiguous bases in the reference
                yield replace_with_n(query, idx_set)
                # Change to list because strings are immutable
                qchars = list(query)
                # Replace the base at each mismatch position with an
                # ambiguous base specifying all possibilities BUT the one
                # we see.
                for idx in idx_set:
                    qchars[idx] = AMBIGUOUS_BASES_COMPLEMENT[qchars[idx]]
                    # Expand to all possibilities for mismatching at this
                    # particular set of positions
                for query_with_mismatches in deambiguate(qchars):
                    yield query_with_mismatches

class PartialMatcher(Matcher):
    def __init__(self, queryset, min_length):
        super().__init__(queryset)
        self.min_length = min_length

        self.partial_queries_left = []
        self.partial_queries_right = []
        for query in self.queryset:
            for s in partial_seqs_left(query, self.min_length):
                self.partial_queries_left.append(s)
            for s in partial_seqs_right(query, self.min_length):
                self.partial_queries_right.append(s)

    def find_match(self, seq):
        for left_partial_query in self.partial_queries_left:
            if seq.startswith(left_partial_query):
               return 0
        for right_partial_query in self.partial_queries_right:
           if seq.endswith(right_partial_query):
               return len(seq) - len(right_partial_query)
        return None

def partial_seqs_left(seq, min_length):
    if min_length < len(seq):
        max_start_idx = len(seq) - min_length + 1
        for start_idx in range(1, max_start_idx):
            yield seq[start_idx:]

def partial_seqs_right(seq, min_length):
    if min_length < len(seq):
        max_remove = len(seq) - min_length
        for num_to_remove in range(1, max_remove + 1):
            yield seq[:-num_to_remove]

def replace_with_n(seq, idxs):
    chars = list(seq)
    for idx in idxs:
        chars[idx] = "N"
    return "".join(chars)

AMBIGUOUS_BASES = {
    "T": "T",
    "C": "C",
    "A": "A",
    "G": "G",
    "R": "AG",
    "Y": "TC",
    "M": "CA",
    "K": "TG",
    "S": "CG",
    "W": "TA",
    "H": "TCA",
    "B": "TCG",
    "V": "CAG",
    "D": "TAG",
    "N": "TCAG",
    }


# Ambiguous base codes for all bases EXCEPT the key
AMBIGUOUS_BASES_COMPLEMENT = {
    "T": "V",
    "C": "D",
    "A": "B",
    "G": "H",
    }


def deambiguate(seq):
    nt_choices = [AMBIGUOUS_BASES[x] for x in seq]
    return ["".join(c) for c in itertools.product(*nt_choices)]


COMPLEMENT_BASES = {
    "T": "A",
    "C": "G",
    "A": "T",
    "G": "C",
    }


def reverse_complement(seq):
    rc = [COMPLEMENT_BASES[x] for x in seq]
    rc.reverse()
    return ''.join(rc)
def filter_paired_main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "input_fastq_fwd", type=argparse.FileType('r'),
        help="Input FASTQ file, forward direction")
    p.add_argument(
        "input_fastq_rev", type=argparse.FileType('r'),
        help="Input FASTQ file, reverse direction")
    p.add_argument(
        "output_fastq_fwd", type=argparse.FileType('w'),
        help="Output FASTQ file, forward direction")
    p.add_argument(
        "output_fastq_rev", type=argparse.FileType('w'),
        help="Output FASTQ file, reverse direction")
    p.add_argument(
        "--min_length", type=int, default=50,
        help="Minimum length for read in either direction")
    args = p.parse_args(argv)

    fwd_reads = FastqRead.parse(args.input_fastq_fwd)
    rev_reads = FastqRead.parse(args.input_fastq_rev)

    for fwd_read, rev_read in zip(fwd_reads, rev_reads):
        fwd_too_short = len(fwd_read.seq) < args.min_length
        rev_too_short = len(rev_read.seq) < args.min_length
        if not (fwd_too_short or rev_too_short):
            args.output_fastq_fwd.write(fwd_read.format_fastq())
            args.output_fastq_rev.write(rev_read.format_fastq())


def main(argv=None):
   p = argparse.ArgumentParser()
   p.add_argument(
      "primer",
      help="Primer sequence to be trimmed")
   p.add_argument(
      "-i", "--input_fastq", type=argparse.FileType('r'),
      help="Input FASTQ file to be trimmed (default: standard input)")
   p.add_argument(
      "-o", "--output_fastq", type=argparse.FileType('w'),
      help="Output fastq after trimming (default: standard output)")
   p.add_argument(
      "--log", type=argparse.FileType('w'),
      help="Log file to record location of primers detected (default: not written)")
   p.add_argument(
      "--num_mismatches", type=int, default=0,
      help="Number of mismatches to primer allowed (default: %(default)s)")
   p.add_argument(
      "--min_length", type=int,
      help="Minimum length for partial primer match (default: full length)")
   p.add_argument(
      "--rev_comp", action='store_true',
      help="Include the reverse complement the primer sequences to trim as well")
   args = p.parse_args(argv)

   if args.input_fastq is None:
      args.input_fastq = sys.stdin

   if args.output_fastq is None:
      args.output_fastq = sys.stdout

   if args.min_length is None:
      args.min_length = len(args.primer)

   if args.rev_comp:
      queryset = deambiguate(args.primer) + deambiguate(reverse_complement(args.primer))
   else:
      queryset = deambiguate(args.primer)

   matchers = [
      CompleteMatcher(queryset, args.num_mismatches),
      PartialMatcher(queryset, args.min_length),
      ]

   for read in FastqRead.parse(args.input_fastq):
       for m in matchers:
           idx = m.find_match(read.seq)
           if idx is not None:
               break
       trimmed_read = read.trim(idx)
       args.output_fastq.write(trimmed_read.format_fastq())

       if args.log:
          if idx is None:
             primer_loc = "NA"
          else:
             primer_loc = str(idx)
          args.log.write("{0}\t{1}\n".format(trimmed_read.desc, primer_loc))
