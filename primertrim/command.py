import argparse
import os
import sys
import tempfile

from .fastq import TrimmableReads
from .matcher import (
    CompleteMatcher, PartialMatcher, AlignmentMatcher,
)
from .dna import deambiguate

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "primer", nargs="+",
        help="Primer sequence to be trimmed")

    io_group = p.add_argument_group("File I/O")
    io_group.add_argument(
        "-i", "--input-fastq", type=argparse.FileType('r'),
        help="Input FASTQ file to be trimmed (default: standard input)")
    io_group.add_argument(
        "-o", "--output-fastq", type=argparse.FileType('w'),
        help="Output fastq after trimming (default: standard output)")
    io_group.add_argument(
        "--log", type=argparse.FileType('w'),
        help="Log file of primers and location (default: not written)")

    complete_group = p.add_argument_group("Complete, partial matching stages")
    complete_group.add_argument(
        "--no-revcomp", action='store_true',
        help=(
            "Don't match the reverse complement during the complete and "
            "partial matching stages"))
    complete_group.add_argument(
        "--mismatches", type=int, default=1,
        help=(
            "Number of mismatches to primer allowed during the complete "
            "matching stage (default: %(default)s)"))
    complete_group.add_argument(
        "--min-partial", type=int, default=8,
        help=(
            "Minimum length of match during the partial matching stage "
            "(default: %(default)s)"))

    alignment_group = p.add_argument_group("Alignment matching stage")
    alignment_group.add_argument(
        "--alignment", action="store_true",
        help="Activate the alignment matching stage")
    alignment_group.add_argument(
        "--alignment-dir",
        help="Directory for alignment files (default: temp directory)")
    alignment_group.add_argument(
        "--threads", type=int,
        help=(
            "Number of CPU threads to use during the alignment stage "
            "(default: all the threads)"))
    args = p.parse_args(argv)

    if args.input_fastq is None:
        args.input_fastq = sys.stdin

    if args.output_fastq is None:
        args.output_fastq = sys.stdout

    queryset = []
    for ambiguous_primer in args.primer:
        for unambiguous_primer in deambiguate(ambiguous_primer):
            queryset.append(unambiguous_primer)

    matchers = [
        CompleteMatcher(queryset, args.mismatches, not args.no_revcomp),
        PartialMatcher(queryset, args.min_partial, not args.no_revcomp),
        ]

    if args.alignment:
        if args.alignment_dir:
            alignment_dir = args.alignment_dir
            if not os.path.exists(alignment_dir):
                os.mkdir(alignment_dir)
        else:
            temp_alignment_dir = tempfile.TemporaryDirectory()
            alignment_dir = temp_alignment_dir.name
        am = AlignmentMatcher(queryset, alignment_dir, args.threads)
        matchers.append(am)

    trimmable_reads = TrimmableReads.from_fastq(args.input_fastq)

    for m in matchers:
        unmatched_seqs = trimmable_reads.get_unmatched_seqs()
        matches_found = m.find_in_seqs(unmatched_seqs)
        for read_id, matchobj in matches_found:
            if matchobj is not None:
                trimmable_reads.register_match(read_id, matchobj)

    for desc, seq, qual in trimmable_reads.get_trimmed_reads():
        if seq != "":
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
