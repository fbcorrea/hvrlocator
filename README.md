# Hypervariable region locator
Workflow to identify spanning hypervariable region(s) from amplicon sequencing variants.
First, the query sequences are aligned to a reference E. coli full-length 16S rRNA gene. That is because the hypervariable region positions are based on this species gene.
Then, the resulting alignment is parsed and the spanning region is idenfied by counting the start and end of the alignment.

## Installing hvreglocator
The easiest way is to install the requirements using conda
```bash
conda create -n hvreglocator mafft biopython
conda activate hvreglocator

git clone https://github.com/felipeborim789/hvreglocator/
```

## Usage
01. Align the reference sequence (ecoli.fa) with the query sequences (query.fa)
```bash
mafft --nuc <(cat ecoli.fa query.fa) > mafft.fa
```

02. Run hvreglocator to infer spanning region
```bash
python3 hvreglocator.py mafft.fa
```
