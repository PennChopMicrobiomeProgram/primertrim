import abc
import collections
import itertools
import os.path

from .dna import (
    AMBIGUOUS_BASES_COMPLEMENT, deambiguate, reverse_complement,
    replace_with_n, partial_seqs_left, partial_seqs_right,
)
from .align import VsearchAligner

PrimerMatch = collections.namedtuple(
    "PrimerMatch", ["method", "start", "mismatches", "primerseq"])


class Matcher(abc.ABC):
    def __init__(self, queryset, match_reverse_complement=True):
        queryset = list(queryset) # We iterate through the queryset twice
        self.queryset = queryset.copy()
        if match_reverse_complement:
            for seq in queryset:
                rc_seq = reverse_complement(seq)
                self.queryset.append(rc_seq)

    def find_in_seqs(self, seqs):
        for seq_id, seq in seqs:
            match = self.find_match(seq)
            yield seq_id, match

    def find_match(self, seq):
        """Returns a PrimerMatch object or None"""
        raise NotImplemented()


class CompleteMatcher(Matcher):
    def __init__(self, queryset, max_mismatch, match_reverse_complement=True):
        super().__init__(queryset, match_reverse_complement)
        self.max_mismatch = max_mismatch

        # To look for near matches in the reference sequence, we
        # generate the set of qeury sequences with 0, 1, or 2 errors
        # and look for these in the reference.  This is our
        # "mismatched queryset."
        possible_mismatches = range(max_mismatch + 1)
        self.mismatched_queryset = [
            self._mismatched_queries(n) for n in possible_mismatches]

    def _mismatched_queries(self, n_mismatch):
        # The generator is provides a sequence for one-time use, but
        # we need to go through it multiple times.  This function
        # wraps the generator to provide a list.
        return list(self._iter_mismatched_queries(n_mismatch))

    def _iter_mismatched_queries(self, n_mismatch):
        # This algorithm is terrible unless the number of mismatches is very small
        assert(n_mismatch in [0, 1, 2, 3])
        for query in self.queryset:
            idx_sets = itertools.combinations(range(len(query)), n_mismatch)
            for idx_set in idx_sets:
                # Replace base at each position with a literal "N", to match
                # ambiguous bases in other sequences
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

    def find_match(self, seq):
        for n_mismatches, queryset in enumerate(self.mismatched_queryset):
            for query in queryset:
                start_idx = seq.find(query)
                if start_idx > -1:
                    end_idx = start_idx + len(query)
                    primerseq = seq[start_idx:end_idx]
                    return PrimerMatch(
                        "Complete", start_idx, n_mismatches, primerseq)


class PartialMatcher(Matcher):
    def __init__(self, queryset, min_length, match_reverse_complement=True):
        super().__init__(queryset, match_reverse_complement)
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
               return PrimerMatch("Partial", 0, 0, left_partial_query)
        for right_partial_query in self.partial_queries_right:
           if seq.endswith(right_partial_query):
               start_idx = len(seq) - len(right_partial_query)
               return PrimerMatch("Partial", start_idx, 0, right_partial_query)


class AlignmentMatcher(Matcher):
    def __init__(self, queryset, alignment_dir, cores=1):
        self.queryset = queryset
        assert(os.path.exists(alignment_dir))
        assert(os.path.isdir(alignment_dir))
        self.alignment_dir = alignment_dir
        self.cores = cores

    def _make_fp(self, filename):
        return os.path.join(self.alignment_dir, filename)

    def find_in_seqs(self, seqs):
        seqs = dict(seqs)
        # Create the file paths
        subject_fp = self._make_fp("subject.fa")
        query_fp = self._make_fp("query.fa")
        result_fp = self._make_fp("vsearch_hits.txt")

        # The database contains the primer sequences
        with open(subject_fp, "w") as f:
            for n, seq in enumerate(self.queryset):
                f.write(">seq{}\n{}\n".format(n, seq))

        a = VsearchAligner(subject_fp)
        hits = a.search(
            seqs.items(), query_fp, result_fp,
            min_id=0.7, threads=self.cores)

        for hit in hits:
            seq_id = hit["qseqid"]
            mismatches = hit["mismatch"] + hit["gapopen"]
            # Vsearch indexes positions starting with 1
            start_idx = hit["qstart"] - 1
            end_idx = hit["qend"]
            assert(start_idx < end_idx)
            seq = seqs[seq_id]
            primerseq = seq[start_idx:end_idx]
            matchobj = PrimerMatch(
                "Alignment", start_idx, mismatches, primerseq)
            yield seq_id, matchobj
