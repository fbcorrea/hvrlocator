# HVRLocator

HVRLocator is a workflow to identify spanning hypervariable region(s) from amplicon sequencing variants or SRA public runs (SRR). It aligns query sequences to a reference E. coli full-length 16S rRNA gene and identifies the spanning region through alignment.

## Using Singularity to run HVRLocator (Recommended)


1. Download the singularity image (*hvrlocator.sif*) to your local folder from the following location:
Go to <https://cloud.sylabs.io/library/jsaraiva/repo/hvrlocator>

**OR**

Paste the following command in your terminal:
```
singularity pull --arch amd64 library://jsaraiva/repo/hvrlocator:hvr
```

## Usage

HVRLocator can process both SRA accession numbers and FASTA files containing ASV sequences. For specific help on each function type one of the following:

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

To use the Random Forest model to also predict the presence of a primer please use the following (also applicable to the processing of a list of runs):

```bash
singularity exec hvrlocator.sif hvrlocator sra -r SRR1585194 --ecoli /path/to/ecoli.fa -o /path/to/output/folder -m /path/to/rf_model.pkl
```

To process a list of SRA runs:

```bash
singularity exec hvrlocator.sif hvrlocator sra -l /path/to/list.txt -o /path/to/output/folder
```
You can specify the location of the E. coli reference file if it's not in the default location:

```bash
singularity exec hvrlocator.sif hvrlocator sra -l /path/to/list.txt --ecoli /path/to/ecoli.fa -o /path/to/output/folder
```

**Note**: The list of SRA runs should be 1 per line.


### Processing ASV FASTA Files

To process a FASTA file containing ASV sequences:

```bash
hvrlocator fasta -f path/to/your/asv_sequences.fasta --ecoli /path/to/ecoli.fa -o /path/to/output/folder
```

### To modify the coverage threshold (default = 0.6) add the "*-t*" flag at the end pf the command (e.g. -t 0.7)" 

## Local Installation (Alternative)

To install HVRLocator locally, follow these steps:

1. Create a new conda environment:

```bash
mamba create -n <ENV_NAME> -y -c bioconda -c conda-forge python=3.9 sra-tools fastp biopython numpy scipy vsearch
```
***Note***: Replace the **ENV_NAME** with a name of your choosing.

2. Activate the environment, clone the repository, and install the package:

```bash
source activate <ENV_NAME> && \
cd <PATH to FOLDER WHERE YOU WANT THE GITHUB REPO TO BE LOCATED> && \
git clone https://github.com/fbcorrea/hvrlocator && \
cd hvrlocator && \
pip install -e .
mamba install -c bioconda -c conda-forge mafft scikit-learn==1.1.3 joblib
```


## Usage

HVRLocator can process both SRA accession numbers and FASTA files containing ASV sequences.

### Processing SRA Accession Numbers

To process an SRA run:

```bash
hvrlocator sra -r SRR1585194 -o /path/to/output_folder
```

You can specify the location of the E. coli reference file if it's not in the default location:

```bash
hvrlocator sra -r SRR1585194 --ecoli /path/to/ecoli.fa -o /path/to/output_folder
```

To use the Random Forest model to also predict the presence of a primer please use the following:

```bash
hvrlocator sra -r SRR1585194 --ecoli /path/to/ecoli.fa -o /path/to/output_folder -m /path/to/rf_model.pkl
```

To process a list of SRA runs (don't forget to add the **-m /path/to/rf_model.pkl** if you wish to also predict primer presence):

```bash
hvrlocator sra -l /path/to/SRA_list.txt -o /path/to/output_folder
```

### Processing ASV FASTA Files

To process a FASTA file containing ASV sequences:

```bash
hvrlocator fasta -f path/to/your/asv_sequences.fasta
```

## Output
The following columns are shown in the output table:

**•	Sample_ID:** Identifier of the processed sample.<br/>
**•	Primer_Presence:** "Yes", "No", or "NA" depending on model output and input quality.<br/>
**•	Score_Primer_Presence:** Probability output from the Random Forest model.<br/>
**•	Median/Avg/Min/Max Alignment Start/End:** Various statistics on read alignment positions.<br/>
**•	Predicted HV Region:** Based on alignment range irrespective of threshold.<br/>
**•	Coverage-based HV Region:** Based on which V-regions passed the specified coverage threshold.<br/>
**•	Coverage_HV_region:** Coverage value of V-regions.<br/>
**•	Warnings:** Alerts about low coverage regions.<br/>
**•	Cov_V1 to Cov_V9:** Coverage values (0-1) for each HV region.<br/>

## Random Forest Model

A detailed description of the Random Forest model generation is available [here](https://github.com/fbcorrea/hvrlocator/blob/main/Random%20Forest%20Model%20for%20Primer%20Presence%20Prediction.docx).

## Project Structure

- `hvrlocator.py`: The main script that handles both SRA and FASTA processing.
- `setup.py`: Used for installing the package.

## Possible Errors and Troubleshooting
**1. FastQ Files Not Found** <br/>
***Error Message:** Error: No FASTQ files found after fastq-dump*<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Ensure that the SRA Run ID is valid.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Check your internet connection for downloading SRA data.<br/>
**2. Low Read Count**<br/>
***Error Message:** Error: Run ID has less than 500 reads and the current sample was skipped*<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Some SRA runs may have low-quality reads or failed sequencing runs.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Consider increasing the read limit in fastq-dump.<br/>
**3. Alignment Issues**<br/>
***Error Message:*** *Error in alignment for ``<ID>``*<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Ensure mafft is installed and available in the system path.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Check if the reference FASTA file (ecoli.fa) is correctly formatted.<br/>
**4. Coverage Too Low for HV Region Assignment**<br/>
***Error Message:** No valid alignment positions*<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Reads may not align properly due to sequencing quality or reference differences.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Check whether trimming parameters in fastp are too strict.<br/>
**5. Missing Columns in TSV Processing**<br/>
***Error Message:** Error: Input TSV file is missing required columns*<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Ensure the input TSV file matches the expected column structure.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Check if the file was manually edited and lost required fields.<br/>
**6. Model Prediction Issues**<br/>
***Error Message:** Model error: No module named ``'sklearn'``*<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Ensure your container or environment includes scikit-learn.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;•	Verify the Random Forest model path with --model is correct and readable. <br/>



## Contributing

Contributions to HVRLocator are welcome. Please feel free to submit a Pull Request.

## License

This project is licensed under the terms of the GNU General Public License v3.0 license.
