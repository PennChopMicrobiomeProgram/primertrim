import itertools

def reverse_complement(seq):
    rc = [COMPLEMENT_BASES[x] for x in seq]
    rc.reverse()
    return ''.join(rc)


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
