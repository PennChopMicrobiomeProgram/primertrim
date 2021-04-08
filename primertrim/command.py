import sys
import argparse
import os

from .fastq import FastqRead
from .matcher import CompleteMatcher, PartialMatcher
from .dna import reverse_complement, deambiguate


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
