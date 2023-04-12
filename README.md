# Primer trim

Detect short primer sequences in FASTQ reads and trim the reads accordingly.

## Installation

```bash
git clone https://github.com/PennChopMicrobiomeProgram/primertrim.git
cd primertrim
conda env create -f primertrim_env.yml -n primertrim_env
conda activate primertrim
pip install .
```

## Algorithm

Primer detection proceeds in three stages: complete matching, partial
matching and (optionally) matching by alignment.

In the complete matching stage, we look for the complete primer
sequence in each read. This stage is implemented in Python, and is
meant to clear out the "easy" matches before the alignment stage. The
user can specify how many mismatches we allow (1 by default). The
algorithm used to detect mismatches starts to get slow for more than 2
mismatches, therefore 3 is the maximum number of mismatches allowed in
this stage.

In the partial matching stage, we try to detect the primer sequence if
it is hanging off the end of the read. The user can specify the
minimum length to signify detection of the partial primer sequence (8
base pairs by default).

For the complete and partial matching stages, the user can specify
whether we search for the reverse complement (yes, by default).

Optionally, the program proceeds to a third stage of matching by
alignment. Here, we use the `vsearch` aligner to detect primer
sequences by semi-global alignment. In this stage, we always search
for the reverse complement. If the user provides an alignment
directory, the alignment files are kept for inspection or
debugging. The `vsearch` program must be installed for the alignment
stage to work.

## Example

```bash
ptrim GCATCGATGAAGAACGCAGC -i sample.fastq -o sample_trimmed.fastq \
    --log sample_trimmed.log --alignment
```
