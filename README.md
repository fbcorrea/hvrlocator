# HVR Locator

HVR Locator is a workflow to identify spanning hypervariable region(s) from amplicon sequencing variants or SRA public runs (SRR). It aligns query sequences to a reference E. coli full-length 16S rRNA gene and identifies the spanning region through alignment.

## Using Singularity to run HVRlocator (Recommended)


1. Download the singularity image (hvrlocator.sif) to your local folder:

## Usage

HVRegLocator can process both SRA accession numbers and FASTA files containing ASV sequences. For specific help on each function type one of the following:

singularity exec hvrlocator.sif hvrlocator sra --help
singularity exec hvrlocator.sif hvrlocator align --help

### Processing SRA Accession Numbers

To process an SRA run:

```bash
singularity exec hvrlocator.sif hvrlocator sra -r SRR1585194 -o /path/to/output/folder
```

You can specify the location of the E. coli reference file if it's not in the default location:

```bash
singularity exec hvrlocator.sif hvrlocator sra -r SRR1585194 --ecoli /path/to/ecoli.fa -o /path/to/output/folder
```

To process a list of SRA runs:

```bash
singularity exec hvrlocator.sif hvrlocator sra -l /path/to/list.txt -o /path/to/output/folder
```
You can specify the location of the E. coli reference file if it's not in the default location:

```bash
singularity exec hvrlocator.sif hvrlocator sra -l /path/to/list.txt --ecoli /path/to/ecoli.fa -o /path/to/output/folder
```

#Note: The list of SRA runs should be 1 per line.


### Processing ASV FASTA Files

To process a FASTA file containing ASV sequences:

```bash
hvreglocator fasta -f path/to/your/asv_sequences.fasta --ecoli /path/to/ecoli.fa -o /path/to/output/folder
```

## Local Installation (Alternative)

To install HVR Locator locally, follow these steps:

1. Create a new conda environment:

```bash
mamba create --prefix /global/apps/hvreglocator/0.2 -y -c bioconda python=3.9 sra-tools mafft fastp biopython numpy scipy vsearch
```

2. Activate the environment, clone the repository, and install the package:

```bash
source activate /global/apps/hvreglocator/0.2 && \
cd /global/apps/hvreglocator/0.2 && \
git clone https://github.com/joao1980/hvrlocator.git && \
cd hvrlocator && \
pip install -e .
```

Note: Replace the GitHub URL with the appropriate URL for your repository.

## Usage

HVRegLocator can process both SRA accession numbers and FASTA files containing ASV sequences.

### Processing SRA Accession Numbers

To process an SRA run:

```bash
hvreglocator sra -r SRR1585194
```

You can specify the location of the E. coli reference file if it's not in the default location:

```bash
hvreglocator sra -r SRR1585194 --ecoli /path/to/ecoli.fa
```

### Processing ASV FASTA Files

To process a FASTA file containing ASV sequences:

```bash
hvreglocator fasta -f path/to/your/asv_sequences.fasta
```

## Output

The script will output the alignment start and end positions, as well as the identified hypervariable region span (e.g., V1-V3).

## Project Structure

- `hvreglocator.py`: The main script that handles both SRA and FASTA processing.
- `setup.py`: Used for installing the package.

## Troubleshooting

If you encounter any issues with finding the E. coli reference file, make sure it's in the same directory as the `hvreglocator.py` script, or use the `--ecoli` argument to specify its location.

## Contributing

Contributions to HVRegLocator are welcome. Please feel free to submit a Pull Request.

## License

This project is licensed under the terms of the MIT license.
