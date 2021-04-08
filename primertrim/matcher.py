import itertools

from .dna import (
    AMBIGUOUS_BASES_COMPLEMENT, deambiguate,
)

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
