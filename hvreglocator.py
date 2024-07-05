import argparse
import statistics
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

def process_sra(input_id, output_dir, ecoli_fa):
    print(f"Processing SRA data for ID: {input_id}")
    working_directory = os.path.join(output_dir, input_id)
    print(f"Working directory: {working_directory}")
    print(f"E. coli reference file: {ecoli_fa}")

    # Ensure working directory exists
    os.makedirs(working_directory, exist_ok=True)

    # Define file paths
    fastq_dump_file_1 = os.path.join(working_directory, f"{input_id}_1.fastq")
    fastq_dump_file_2 = os.path.join(working_directory, f"{input_id}_2.fastq")
    fastp_file_1 = os.path.join(working_directory, f"{input_id}_fastp_1.fastq")
    fastp_file_2 = os.path.join(working_directory, f"{input_id}_fastp_2.fastq")
    merged_file = os.path.join(working_directory, f"{input_id}_merged.fasta")
    sed_file = os.path.join(working_directory, f"{input_id}_processed.fasta")
    aligned_file = os.path.join(working_directory, f"{input_id}_aligned.fasta")
    final_output_file = os.path.join(working_directory, f"{input_id}_hvreglocator_output.txt")

    print(f"Output files will be saved in: {working_directory}")

    # Step 1: Run fastq-dump to get reads (first 1000)
    print("\nStep 1: Running fastq-dump...")
    fastq_dump_command = f"fastq-dump  --split-files -X 1000 {input_id} -O {working_directory}"
    if not run_command(fastq_dump_command):
        print("Error in fastq-dump command. Exiting.")
        return False

    # Debug: List files after fastq-dump
    print("\nListing files after fastq-dump:")
    run_command(f"ls -l {working_directory}")

    # Check if we have paired-end or single-end reads
    is_paired = os.path.exists(fastq_dump_file_2)
    print(f"Is paired-end: {is_paired}")
    print(f"File 1 exists: {os.path.exists(fastq_dump_file_1)}")
    print(f"File 2 exists: {os.path.exists(fastq_dump_file_2)}")

    # Step 2: Run fastp
    print("\nStep 2: Running fastp...")
    if is_paired:
        print("Processing paired-end reads with fastp...")
        fastp_command = f"fastp --in1 {fastq_dump_file_1} --in2 {fastq_dump_file_2} --out1 {fastp_file_1} --out2 {fastp_file_2} --detect_adapter_for_pe"
    else:
        print("Processing single-end reads with fastp...")
        fastp_command = f"fastp --in1 {fastq_dump_file_1} --out1 {fastp_file_1}"

    if not run_command(fastp_command):
        print("Error in fastp command. Exiting.")
        return False

    # Debug: List files after fastp
    print("\nListing files after fastp:")
    run_command(f"ls -l {working_directory}")

    # Step 3: Run vsearch for paired-end reads or sed for single-end reads
    if is_paired:
        print("\nStep 3: Running vsearch to merge paired-end reads...")
        vsearch_command = f"vsearch --fastq_mergepairs {fastp_file_1} --reverse {fastp_file_2} --fastaout {merged_file}"
        if not run_command(vsearch_command):
            print("Error in vsearch command. Exiting.")
            return False
        sed_file = merged_file
    else:
        print("\nStep 3: Running sed to convert FASTQ to FASTA...")
        sed_command = f"sed -n '1~4s/^@/>/p;2~4p' {fastp_file_1} > {sed_file}"
        if not run_command(sed_command):
            print("Error in sed command. Exiting.")
            return False

    # Debug: List files after vsearch/sed
    print("\nListing files after vsearch/sed:")
    run_command(f"ls -l {working_directory}")

    # Step 4: Align reads to reference
    print("\nStep 4: Aligning reads to reference...")
    align_command = f"mafft --thread 8 --6merpair --nuc --keeplength --addfragments {sed_file} {ecoli_fa} > {aligned_file}"
    if not run_command(align_command):
        print("Error in alignment. Exiting.")
        return False

    # Step 5: Process aligned file
    print("\nStep 5: Processing aligned file...")
    hvreglocator_output = process_fasta(aligned_file)

    with open(final_output_file, 'w') as f:
        f.write(hvreglocator_output)

    print(f"\nProcess completed successfully. Final output saved in: {final_output_file}")
    return final_output_file

def run_command(command):
    print(f"Executing command: {command}")
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print("Command executed successfully.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command. Return code: {e.returncode}")
        print(f"Error message: {e.stderr}")
        return False
    except Exception as e:
        print(f"Exception occurred: {str(e)}")
        return False

def process_fasta(input_fasta):
    print(f"Processing FASTA file: {input_fasta}")
    fasta_sequences = list(SeqIO.parse(input_fasta, 'fasta'))
    print(f"Total sequences read: {len(fasta_sequences)}")
    
    if len(fasta_sequences) < 2:
        print("Error: Not enough sequences in the alignment file.")
        return "Error\tInsufficient sequences for analysis"
    
    start_positions = []
    end_positions = []
    
    for seq in fasta_sequences:
        seq_str = str(seq.seq)
        # Find start position (first non-gap character)
        start_pos = next((i for i, char in enumerate(seq_str) if char != '-'), None)
        # Find end position (last non-gap character)
        end_pos = len(seq_str) - next((i for i, char in enumerate(reversed(seq_str)) if char != '-'), None) - 1
        
        if start_pos is not None and end_pos is not None:
            start_positions.append(start_pos)
            end_positions.append(end_pos)
    
    if not start_positions or not end_positions:
        print("Error: No valid alignment positions found.")
        return "Error\tNo valid alignment positions"
    
    # Calculate median start and end positions
    median_start = int(statistics.median(start_positions))
    median_end = int(statistics.median(end_positions))
    
    # Determine hypervariable regions
    region_start = determine_region(median_start)
    region_end = determine_region(median_end, is_start=False)
    
    # Create TSV output
    headers = ["Median_Alignment_start", "Median_Alignment_end", "HV_region_start", "HV_region_end"]
    values = [median_start, median_end, f"V{region_start}", f"V{region_end}"]
    results = "\t".join(headers) + "\n" + "\t".join(map(str, values))
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

def process_id_list(id_list_file, output_dir, ecoli_fa):
    with open(id_list_file, 'r') as f:
        run_ids = f.read().splitlines()
    
    for run_id in run_ids:
        print(f"\nProcessing run ID: {run_id}")
        process_sra(run_id, output_dir, ecoli_fa)

def main():
    print("Starting HVRegLocator...")
    parser = argparse.ArgumentParser(description="HVRegLocator: Process SRA data or FASTA file to locate hypervariable regions.")
    parser.add_argument("command", choices=["sra", "fasta"], help="Choose 'sra' to process SRA data or 'fasta' to process a FASTA file")
    parser.add_argument("-r", "--run-id", help="SRA Run ID for processing (required for 'sra' command unless -l is used)")
    parser.add_argument("-l", "--list-file", help="File containing list of SRA Run IDs to process (one per line)")
    parser.add_argument("-f", "--fasta-file", help="Input FASTA file (required for 'fasta' command)")
    parser.add_argument("-e", "--ecoli", help="Path to ecoli.fa file", default=None)
    parser.add_argument("-o", "--output-dir", help="Output directory", required=True)
    
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

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    if args.command == "sra":
        if args.list_file:
            print(f"Processing list of SRA run IDs from file: {args.list_file}")
            process_id_list(args.list_file, args.output_dir, ecoli_fa)
        elif args.run_id:
            print(f"Processing single SRA run ID: {args.run_id}")
            process_sra(args.run_id, args.output_dir, ecoli_fa)
        else:
            print("Error: The 'sra' command requires either --run-id or --list-file argument")
            parser.error("The 'sra' command requires either --run-id or --list-file argument")
    
    elif args.command == "fasta":
        if not args.fasta_file:
            print("Error: The 'fasta' command requires the --fasta-file argument")
            parser.error("The 'fasta' command requires the --fasta-file argument")
        
        print(f"Processing FASTA file: {args.fasta_file}")
        process_fasta(args.fasta_file)

    print("HVRegLocator execution complete.")

if __name__ == "__main__":
    main()