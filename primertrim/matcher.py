import abc
import collections
import itertools

from .dna import (
    AMBIGUOUS_BASES_COMPLEMENT, deambiguate, replace_with_n,
    partial_seqs_left, partial_seqs_right,
)


PrimerMatch = collections.namedtuple(
    "PrimerMatch", ["method", "start", "mismatches", "primerseq"])


class Matcher(abc.ABC):
    def __init__(self, queryset):
        self.queryset = queryset

    def find_in_seqs(self, seqs):
        for seq_id, seq in seqs:
            match = self.find_match(seq)
            yield seq_id, match

    def find_match(self, seq):
        """Returns a PrimerMatch object or None"""
        raise NotImplemented()


class CompleteMatcher(Matcher):
    def __init__(self, queryset, max_mismatch):
        super().__init__(queryset)
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
               return PrimerMatch("Partial", 0, 0, left_partial_query)
        for right_partial_query in self.partial_queries_right:
           if seq.endswith(right_partial_query):
               start_idx = len(seq) - len(right_partial_query)
               return PrimerMatch("Partial", start_idx, 0, right_partial_query)
