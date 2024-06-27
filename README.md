# Hypervariable region locator

Workflow to identify spanning hypervariable region(s) from amplicon sequencing variants
or SRA public runs (SRR).
First, the query sequences are aligned to a reference E. coli full-length 16S rRNA gene. 
That is because the hypervariable region positions are based on this species gene.
Then, the resulting alignment the spanning region is idenfied through alignment.

# Installing hvreglocator

The easiest way is to install the requirements using Conda (in this case through Mamba)
```bash
mamba create -n hvreglocator -c bioconda sra-tools mafft biopython fastp
source activate hvreglocator
git clone https://github.com/felipeborim789/hvreglocator/
```

# Scripts

- **hvreglocator.py** works for resolved ASVs. The input argument expects the path to a fasta file containing a single ASV DNA sequence.

- **hvreglocator_mod.py** module that works for SRA accession numbers. The input argument expects a public SRR accession number (SRA Run containing the sequencing data).

- **main.py** calls the whole pipeline for SRA and runs **hvreglocator_mod.py**.

## Usage with SRA accesion numbers

Example with run accession SRR1585194.

According to the name of the run this is V1-V3:
Illumina MiSeq paired end sequencing; M0-inoculate Illumina 16s V1-V3 Environmental Sample


```bash
python3 main.py -i SRR1585194
```

**Results**
Altough the expected output is V1-V3, I am getting V1-V2. I need to check why.


## Usage with ASVs

01. Align the reference sequence (ecoli.fa) with the query sequences (query.fa)
```bash
mafft --nuc --keeplength --addfragments query.fa ecoli.fa > mafft.fa
```
Info about "--addfragments" and "--keeplength" options: https://mafft.cbrc.jp/alignment/server/add.html

02. Run hvreglocator to infer spanning region
```bash
python3 hvreglocator.py mafft.fa
```
