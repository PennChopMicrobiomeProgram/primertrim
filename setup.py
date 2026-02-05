from setuptools import setup

setup(
    name="primertrim",
    version="0.0.3",
    description="Trim primer sequences from FASTQ files",
    author="PennCHOP Microbiome Program",
    author_email="BITTINGERK@chop.edu",
    url="https://github.com/PennChopMicrobiomeProgram/primertrim",
    packages=["primertrim"],
    entry_points={
        "console_scripts": [
            "ptrim=primertrim.command:main",
        ],
    },
)
