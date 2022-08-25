# Hypervariable region locator
Pipeline to identify spanning hypervariable region(s) from amplicon sequencing variants.

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
