from setuptools import setup

setup(
    name='primertrim',
    version='0.0.1',
    description='Trim primer sequences from FASTQ files',
    author='PennCHOP Microbiome Program',
    author_email='BITTINGERK@chop.edu',
    url='https://github.com/PennChopMicrobiomeProgram/Primer_trim',
    packages=['primertrim'],
    entry_points = {
        'console_scripts': [
            'remove_primers.py=primertrim.command:main',
        ],
    }
)
