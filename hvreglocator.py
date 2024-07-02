import argparse
import subprocess
import os
import sys
from Bio import SeqIO

def run_command(command, working_directory=None):
    print(f"Executing command: {command}")
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=working_directory)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Error executing command. Return code: {process.returncode}")
            print(f"Error message: {stderr.decode('utf-8')}")
            return False
        print("Command executed successfully.")
        return True
    except Exception as e:
        print(f"Exception occurred: {str(e)}")
        return False

def process_sra(input_id, working_directory, ecoli_fa):
    print(f"Processing SRA data for ID: {input_id}")
    print(f"Working directory: {working_directory}")
    print(f"E. coli reference file: {ecoli_fa}")

    # Ensure working directory exists
    os.makedirs(working_directory, exist_ok=True)

    # Change to the working directory
    os.chdir(working_directory)

    # Define file paths (now relative to working_directory)
    fastq_dump_file_1 = f"{input_id}_1.fastq"
    fastq_dump_file_2 = f"{input_id}_2.fastq"
    merged_file = f"{input_id}_merged.fasta"
    aligned_file = f"{input_id}_aligned.fasta"
    final_output_file = "hvreglocator_output.txt"

    print(f"Output files will be saved in: {working_directory}")

    # Step 1: Run fastq-dump to get paired-end reads
    print("\nStep 1: Running fastq-dump...")
    fastq_dump_command = f"fastq-dump --split-files {input_id}"
    if not run_command(fastq_dump_command, working_directory):
        print("Error in fastq-dump command. Exiting.")
        return False

    # Step 2: Merge paired-end reads using vsearch and output as FASTA
    print("\nStep 2: Merging paired-end reads...")
    merge_command = f"vsearch --fastq_mergepairs {fastq_dump_file_1} --reverse {fastq_dump_file_2} --fastaout {merged_file}"
    if not run_command(merge_command, working_directory):
        print("Error in merging paired-end reads. Exiting.")
        return False

    # Step 3: Align merged reads to reference
    print("\nStep 3: Aligning merged reads to reference...")
    align_command = f"mafft --thread 8 --6merpair --nuc --keeplength --addfragments {merged_file} {ecoli_fa} > {aligned_file}"
    if not run_command(align_command, working_directory):
        print("Error in alignment. Exiting.")
        return False

    # Step 4: Process aligned file
    print("\nStep 4: Processing aligned file...")
    hvreglocator_output = process_fasta(aligned_file)
    
    with open(final_output_file, 'w') as f:
        f.write(hvreglocator_output)

    print(f"\nProcess completed successfully. Final output saved in: {os.path.join(working_directory, final_output_file)}")
    return os.path.join(working_directory, final_output_file)

def process_fasta(input_fasta):
    print(f"Processing FASTA file: {input_fasta}")
    
    fasta_sequences = list(SeqIO.parse(input_fasta, 'fasta'))
    print(f"Total sequences read: {len(fasta_sequences)}")
    
    if len(fasta_sequences) < 2:
        print("Error: Not enough sequences in the alignment file.")
        return "Error: Insufficient sequences for analysis."

    ref_seq = fasta_sequences[0]
    query_seq = fasta_sequences[1]

    ref_str = str(ref_seq.seq)
    query_str = str(query_seq.seq)

    # Find start position (first non-gap character in query)
    start_pos = next((i for i, char in enumerate(query_str) if char != '-'), None)

    # Find end position (last non-gap character in query)
    end_pos = len(query_str) - next((i for i, char in enumerate(reversed(query_str)) if char != '-'), None) - 1

    # Determine hypervariable regions
    region_start = determine_region(start_pos)
    region_end = determine_region(end_pos, is_start=False)

    results = f"\nResults:\n"
    results += f"Alignment start: {start_pos}\n"
    results += f"Alignment end: {end_pos}\n"
    results += f"Hypervariable region span: V{region_start}-V{region_end}\n"
    
    print(results)
    return results

def determine_region(pos, is_start=True):
    regions = [
        (1, 0, 69),
        (2, 99, 157),
        (3, 227, 440),
        (4, 500, 590),
        (5, 650, 828),
        (6, 857, 1000),
        (7, 1036, 1119),
        (8, 1158, 1243),
        (9, 1295, 1435)
    ]

    if not is_start:
        regions.reverse()

    for region, start, end in regions:
        if start <= pos <= end:
            return region

    return "Unknown"

def main():
    print("Starting HVRegLocator...")
    parser = argparse.ArgumentParser(description="HVRegLocator: Process SRA data or FASTA file to locate hypervariable regions.")
    parser.add_argument("command", choices=["sra", "fasta"], help="Choose 'sra' to process SRA data or 'fasta' to process a FASTA file")
    parser.add_argument("-r", "--run-id", help="SRA Run ID for processing (required for 'sra' command)")
    parser.add_argument("-f", "--fasta-file", help="Input FASTA file (required for 'fasta' command)")
    parser.add_argument("-e", "--ecoli", help="Path to ecoli.fa file", default=None)
    parser.add_argument("-w", "--working-dir", help="Working directory", default=os.getcwd())
    
    args = parser.parse_args()
    print(f"Command-line arguments: {args}")

    # Determine the path to ecoli.fa
    if args.ecoli:
        ecoli_fa = os.path.abspath(args.ecoli)
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
        working_directory = os.path.abspath(args.working_dir)
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