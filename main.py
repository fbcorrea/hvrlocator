import argparse
import subprocess
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process the data pipeline with a given input ID.")
parser.add_argument("-i", "--input_id", required=True, help="Input ID for the data processing pipeline.")
args = parser.parse_args()

# Create Wdir if it doesn't exist
working_directory = "Wdir"
if not os.path.exists(working_directory):
    os.makedirs(working_directory)

def run_command(command, input_file=None, output_file=None):
    infile = open(input_file, 'r') if input_file else None
    outfile = open(output_file, 'w') if output_file else None

    try:
        process = subprocess.Popen(command, stdin=infile, stdout=outfile if outfile else subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        _, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Error: {stderr.decode('utf-8')}")
            return False
    finally:
        if infile:
            infile.close()
        if outfile:
            outfile.close()

    return True

# Define file paths with Wdir directory
fastq_dump_file = os.path.join(working_directory, "fastq_dump_output.txt")
fastp_file = os.path.join(working_directory, "fastp_output.txt")
sed_file = os.path.join(working_directory, "sed_output.txt")
mafft_file = os.path.join(working_directory, "mafft_output.txt")
final_output_file = os.path.join(working_directory, "hvreglocator_output.txt")

# Define and run each command
print("Running fastq-dump...")
fastq_dump_command = f"fastq-dump --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' -Z {args.input_id} | head -n 4096 > {fastq_dump_file}"
if not run_command(fastq_dump_command):
    print("Error in fastq-dump command")
    exit(1)

print("Running fastp...")
fastp_command = f"cat {fastq_dump_file} | fastp --stdin --stdout --interleaved_in --detect_adapter_for_pe > {fastp_file}"
if not run_command(fastp_command):
    print("Error in fastp command")
    exit(1)

print("Running sed...")
sed_command = f"cat {fastp_file} | sed -n '1~4s/^@/>/p;2~4p' > {sed_file}"
if not run_command(sed_command):
    print("Error in sed command")
    exit(1)

print("Running mafft...")
mafft_command = f"mafft --nuc --keeplength --addfragments {sed_file} ecoli.fa > {mafft_file}"
if not run_command(mafft_command):
    print("Error in mafft command")
    exit(1)

print("Running hvreglocator_mod.py...")
hvreglocator_command = f"python hvreglocator_mod.py {mafft_file} > {final_output_file}"
if not run_command(hvreglocator_command):
    print("Error in hvreglocator_mod.py command")
    exit(1)

print("Process completed successfully. Output saved in:", final_output_file)

