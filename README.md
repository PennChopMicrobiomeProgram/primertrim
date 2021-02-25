# Primer trim

## Installation:
```bash
git clone https://github.com/PennChopMicrobiomeProgram/Primer_trim.git
cd Primer_trim
conda activate (env_name)
pip install -e ./
```

## Arguments:
-**--input_fastq**: Either forward (R1) or reverse (R2) fastq file
-**--output_fastq**: Fastq file (R1/R2) with primer sequence trimmed from input fastq file
-**--log**: A log file that records the position of the trimmed primer sequence
-**--num_mismatches**: Number of mismatches to primer allowed (default: 0)
-**--min_length**: Minimum length for partial primer match (default: full length of primer sequence)
-**--rev_comp**: Include the reverse complement the primer sequences to trim as well

## Example:
```
remove_primers.py ${PRIMER_SEQUENCE} -i ${INPUT_R1_FASTQ_PATH} -o ${OUTPUT_R1_FASTQ_PATH} --log ${LOG_PATH} --num_mismatches ${INT} --min_length ${INT} --rev_comp
```

Remark: both `${INPUT_R1_FASTQ_PATH}` and `${OUTPUT_R1_FASTQ_PATH}` __must__ end with "R1.fastq"