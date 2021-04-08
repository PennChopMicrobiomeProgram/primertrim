import sys
import argparse

from .fastq import FastqRead
from .matcher import CompleteMatcher, PartialMatcher
from .dna import reverse_complement, deambiguate


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
