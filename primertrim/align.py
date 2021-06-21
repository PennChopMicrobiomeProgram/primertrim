import os.path
import subprocess
import tempfile

DEFAULT_BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "qlen", "slen", "qseq", "sseq", "qstrand",
]

BLAST_FIELD_TYPES = {
    "qseqid": str,
    "sseqid": str,
    "pident": float,
    "length": int,
    "mismatch": int,
    "gapopen": int,
    "qstart": int,
    "qend": int,
    "sstart": int,
    "send": int,
    "qlen": int,
    "slen": int,
    "qseq": str,
    "sseq": str,
    "qstrand": str,
}

BLAST_TO_VSEARCH = {
    "qseqid": "query",
    "sseqid": "target",
    "pident": "id2",
    "length": "alnlen",
    "mismatch": "mism",
    "gapopen": "gaps",
    "qstart": "qilo",
    "qend": "qihi",
    "sstart": "tilo",
    "send": "tihi",
    "qlen": "qs",
    "slen": "ts",
    "qseq": "qrow",
    "sseq": "trow",
    "qstrand": "qstrand",
}

def write_fasta(f, seqs):
    for seq_id, seq in seqs:
        f.write(">{}\n{}\n".format(seq_id, seq))

class VsearchAligner:
    def __init__(self, ref_seqs_fp):
        self.ref_seqs_fp = ref_seqs_fp
        self.fields = DEFAULT_BLAST_FIELDS
        self.convert_types = True
        self.stderr = subprocess.DEVNULL

    def search(self, seqs, input_fp=None, output_fp=None, **kwargs):
        """Search seqs and return hits"""
        if input_fp is None:
            infile = tempfile.NamedTemporaryFile(
                suffix=".fasta", mode="w+t", encoding="utf-8")
            write_fasta(infile, seqs)
            infile.seek(0)
            input_fp = infile.name
        else:
            with open(input_fp, "w") as f:
                write_fasta(f, seqs)

        if output_fp is None:
            outfile = tempfile.NamedTemporaryFile(
                suffix=".txt", mode="wt")
            output_fp = outfile.name

        self._call(input_fp, output_fp, **kwargs)

        with open(output_fp) as f:
            for hit in self.parse(f):
                yield hit

    def _call(self, query_fp, output_fp, min_id=0.85, threads=None):
        id_arg = "{:.3f}".format(min_id)
        userfields_arg = "+".join(BLAST_TO_VSEARCH[f] for f in self.fields)
        args = [
            "vsearch",
            "--usearch_global", query_fp,
            "--minseqlength", "10",
            "--mincols", "10",
            "--id", id_arg,
            "--wordlength", "4",
            "--strand", "both",
            "--maxaccepts", "4",
            "--minwordmatches", "3",
            "--top_hits_only",
            "--userfields", userfields_arg,
            "--db", self.ref_seqs_fp,
            "--userout", output_fp,
        ]
        if threads is not None:
            threads_arg = "{:d}".format(threads)
            args.extend(["--threads", threads_arg])
        subprocess.check_call(args, stderr=self.stderr)

    def parse(self, f):
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            res = dict(zip(self.fields, vals))
            if self.convert_types:
                for field in self.fields:
                    fcn = BLAST_FIELD_TYPES[field]
                    res[field] = fcn(res[field])
            if res["qstrand"] == "-":
                qstart_temp = res["qlen"]-res["qend"]+1
                res["qend"] = res["qlen"]-res["qstart"]+1
                res["qstart"] = qstart_temp
            yield res
