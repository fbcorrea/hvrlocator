import argparse
import subprocess
import os
import sys
import numpy as np
from scipy import stats
from Bio import SeqIO
import matplotlib.pyplot as plt

def run_command(command, input_file=None, output_file=None):
    print(f"Executing command: {command}")
    infile = open(input_file, 'r') if input_file else None
    outfile = open(output_file, 'w') if output_file else None
    try:
        process = subprocess.Popen(command, stdin=infile, stdout=outfile if outfile else subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        _, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Error executing command. Return code: {process.returncode}")
            print(f"Error message: {stderr.decode('utf-8')}")
            return False
        print("Command executed successfully.")
    finally:
        if infile:
            infile.close()
        if outfile:
            outfile.close()
    return True

def process_sra(input_id, working_directory, ecoli_fa):
    print(f"Processing SRA data for ID: {input_id}")
    print(f"Working directory: {working_directory}")
    print(f"E. coli reference file: {ecoli_fa}")

    # Define file paths with Wdir directory
    fastq_dump_file = os.path.join(working_directory, "fastq_dump_output.txt")
    fastp_file = os.path.join(working_directory, "fastp_output.txt")
    sed_file = os.path.join(working_directory, "sed_output.txt")
    mafft_file = os.path.join(working_directory, "mafft_output.txt")
    final_output_file = os.path.join(working_directory, "hvreglocator_output.txt")

    print(f"Output files will be saved in: {working_directory}")

    # Define and run each command
    print("\nStep 1: Running fastq-dump...")
    fastq_dump_command = f"fastq-dump --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' -Z {input_id} | head -n 8192 > {fastq_dump_file}"
    if not run_command(fastq_dump_command):
        print("Error in fastq-dump command. Exiting.")
        return False

    print("\nStep 2: Running fastp...")
    fastp_command = f"cat {fastq_dump_file} | fastp --stdin --stdout --interleaved_in --detect_adapter_for_pe > {fastp_file}"
    if not run_command(fastp_command):
        print("Error in fastp command. Exiting.")
        return False

    print("\nStep 3: Running sed...")
    sed_command = f"cat {fastp_file} | sed -n '1~4s/^@/>/p;2~4p' > {sed_file}"
    if not run_command(sed_command):
        print("Error in sed command. Exiting.")
        return False

    print("\nStep 4: Running mafft...")
    mafft_command = f"mafft --thread 8 --6merpair --nuc --keeplength --addfragments {sed_file} {ecoli_fa} > {mafft_file}"
    if not run_command(mafft_command):
        print("Error in mafft command. Exiting.")
        return False

    print("\nStep 5: Running hvreglocator...")
    hvreglocator_output = process_fasta(mafft_file)
    
    with open(final_output_file, 'w') as f:
        f.write(hvreglocator_output)

    print(f"\nProcess completed successfully. Final output saved in: {final_output_file}")
    return final_output_file

def process_fasta(input_fasta):
    print(f"Processing FASTA file: {input_fasta}")
    
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Total sequences read: {len(sequences)}")
    
    ref_seq = sequences.pop(0)
    print(f"Sequences after removing reference: {len(sequences)}")

    filtered_sequences = [record for record in sequences if '/' in record.id]
    print(f"Filtered sequences: {len(filtered_sequences)}")

    start_positions = []
    end_positions = []

    for record in filtered_sequences:
        seq_str = str(record.seq)
        start = next((i for i, char in enumerate(seq_str) if char != '-'), None)
        end = len(seq_str) - next((i for i, char in enumerate(reversed(seq_str)) if char != '-'), None)
        if start is not None and end is not None:
            start_positions.append(start)
            end_positions.append(end)

    start_positions = np.array(start_positions)
    end_positions = np.array(end_positions)
    
    print(f"Start positions: min={np.min(start_positions)}, max={np.max(start_positions)}, mean={np.mean(start_positions):.2f}, median={np.median(start_positions)}")
    print(f"End positions: min={np.min(end_positions)}, max={np.max(end_positions)}, mean={np.mean(end_positions):.2f}, median={np.median(end_positions)}")
    
    # Create histograms
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.hist(start_positions, bins=50)
    plt.title("Start Positions Histogram")
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    
    plt.subplot(1, 2, 2)
    plt.hist(end_positions, bins=50)
    plt.title("End Positions Histogram")
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    
    plt.tight_layout()
    plt.savefig("position_histograms.png")
    print("Histograms saved as 'position_histograms.png'")

    # Use mode (most common value) for start and end
    start_pos = int(np.bincount(start_positions).argmax())
    end_pos = int(np.bincount(end_positions).argmax())

    print(f"Using start_pos: {start_pos}, end_pos: {end_pos}")

    # Adjust these ranges based on the E. coli 16S rRNA gene positions for 515F-Y and 926R
    if 500 <= start_pos < 550:
        region_start = 4  # 515F-Y primer
    else:
        region_start = "Unknown"

    if 900 <= end_pos < 950:
        region_end = 5  # 926R primer
    elif 750 <= end_pos < 900:
        region_end = 4  # In case the alignment ends within V4
    else:
        region_end = "Unknown"

    print(f"Determined start region: V{region_start}")
    print(f"Determined end region: V{region_end}")

    results = f"\nResults:\n"
    results += f"Alignment start: {start_pos}\n"
    results += f"Alignment end: {end_pos}\n"
    results += f"Hypervariable region span: V{region_start}-V{region_end}\n"
    
    print(results)
    return results

def main():
    print("Starting HVRegLocator...")
    parser = argparse.ArgumentParser(description="HVRegLocator: Process SRA data or FASTA file to locate hypervariable regions.")
    parser.add_argument("command", choices=["sra", "fasta"], help="Choose 'sra' to process SRA data or 'fasta' to process a FASTA file")
    parser.add_argument("-r", "--run-id", help="SRA Run ID for processing (required for 'sra' command)")
    parser.add_argument("-f", "--fasta-file", help="Input FASTA file (required for 'fasta' command)")
    parser.add_argument("-e", "--ecoli", help="Path to ecoli.fa file", default=None)
    
    args = parser.parse_args()
    print(f"Command-line arguments: {args}")

    # Determine the path to ecoli.fa
    if args.ecoli:
        ecoli_fa = args.ecoli
    else:
        # Use the directory of the script to find ecoli.fa
        script_dir = os.path.dirname(os.path.abspath(__file__))
        ecoli_fa = os.path.join(script_dir, "ecoli.fa")
    
    print(f"Using E. coli reference file: {ecoli_fa}")

    if not os.path.exists(ecoli_fa):
        print(f"Error: E. coli reference file not found at {ecoli_fa}")
        sys.exit(1)

    if args.command == "sra":
        if not args.run_id:
            print("Error: The 'sra' command requires the --run-id argument")
            parser.error("The 'sra' command requires the --run-id argument")
        
        print(f"Processing SRA data for run ID: {args.run_id}")
        working_directory = "Wdir"
        if not os.path.exists(working_directory):
            print(f"Creating working directory: {working_directory}")
            os.makedirs(working_directory)
        
        output_file = process_sra(args.run_id, working_directory, ecoli_fa)
        if not output_file:
            print("SRA processing failed. Exiting.")
    
    elif args.command == "fasta":
        if not args.fasta_file:
            print("Error: The 'fasta' command requires the --fasta-file argument")
            parser.error("The 'fasta' command requires the --fasta-file argument")
        
        print(f"Processing FASTA file: {args.fasta_file}")
        process_fasta(args.fasta_file)

    print("HVRegLocator execution complete.")

if __name__ == "__main__":
    main()