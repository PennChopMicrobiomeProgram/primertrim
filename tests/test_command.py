from pathlib import Path

from primertrim.command import main


DATA_DIR = Path(__file__).parent / "data"


def data_fp(filename):
    return str(DATA_DIR / filename)


def read_from(filepath):
    with open(filepath) as f:
        res = f.readlines()
    return res

def test_main_script(tmp_path):
    input_fp = data_fp("Sub10003.V1.sputum.redo_R1.fastq")
    output_fp = str(tmp_path / "out.fastq")
    log_fp = str(tmp_path / "out.log")
    args = [
        "GCATCGATGAAGAACGCAGC",
        "-i", input_fp,
        "-o", output_fp,
        "--log", log_fp,
        "--mismatches", "0",
        "--min-partial", "100",
    ]
    main(args)

    expected_log_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_R1.log")
    assert read_from(log_fp) == read_from(expected_log_fp)

    expected_output_fp = data_fp("no_primer_Sub10003.V1.sputum.redo_R1.fastq")
    assert read_from(output_fp) == read_from(expected_output_fp)
