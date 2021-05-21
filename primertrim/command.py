import argparse
import os
import sys
import tempfile

from .trimmable_reads import TrimmableReads
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
        help="Output FASTQ file after trimming (default: standard output)")
    io_group.add_argument(
        "--log", type=argparse.FileType('w'),
        help="Log file of primers and location (default: not written)")
    io_group.add_argument(
        "--min-length", type=int, default=50,
        help=(
            "Minimum length of reads written to the output FASTQ file. "
            "(default: %(default)s)"))

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

    output_reads = trimmable_reads.output_reads(args.min_length)
    write_fastq(args.output_fastq, output_reads)

    output_loginfo = trimmable_reads.output_loginfo()
    write_log(args.log, output_loginfo, trimmable_reads.loginfo_colnames)


def write_fastq(f, reads):
    for desc, seq, qual in reads:
        f.write("@{0}\n{1}\n+\n{2}\n".format(desc, seq, qual))


def write_log(f, loginfo, colnames=None):
    if colnames:
        f.write("\t".join(colnames))
        f.write("\n")
    for vals in loginfo:
        f.write("\t".join(str(v) if v is not None else "" for v in vals))
        f.write("\n")
