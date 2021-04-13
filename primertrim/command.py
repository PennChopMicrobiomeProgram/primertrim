import sys
import argparse
import tempfile

from .fastq import TrimmableReads, parse_fastq
from .matcher import (
    CompleteMatcher, PartialMatcher, AlignmentMatcher,
)
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
    p.add_argument(
        "--skip-alignment", action="store_true",
        help="Skip the alignment stage")
    p.add_argument(
        "--alignment-dir",
        help="Directory for alignment files (default: temp directory)")
    p.add_argument(
        "--threads", type=int,
        help="Number of CPU threads to use during alignment stage")
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

    if not args.skip_alignment:
        if args.alignment_dir:
            alignment_dir = args.alignment_dir
        else:
            temp_alignment_dir = tempfile.TemporaryDirectory()
            alignment_dir = temp_alignment_dir.name
        am = AlignmentMatcher(queryset, alignment_dir, args.threads)
        matchers.append(am)

    trimmable_reads = TrimmableReads(parse_fastq(args.input_fastq))

    for m in matchers:
        unmatched_seqs = trimmable_reads.get_unmatched_seqs()
        matches_found = m.find_in_seqs(unmatched_seqs)
        for read_id, matchobj in matches_found:
            if matchobj is not None:
                trimmable_reads.register_match(read_id, matchobj)

    for desc, seq, qual in trimmable_reads.get_trimmed_reads():
        args.output_fastq.write("@{0}\n{1}\n+\n{2}\n".format(desc, seq, qual))

    if args.log:
        for read_id, matchobj in trimmable_reads.matches.items():
            if matchobj is None:
                args.log.write("{0}\t\t\t\t\n".format(read_id))
            else:
                args.log.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                    read_id, matchobj.method, matchobj.start,
                    matchobj.mismatches, matchobj.primerseq,
                ))
