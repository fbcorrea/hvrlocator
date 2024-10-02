import argparse
import statistics
import subprocess
import os
import sys
from Bio import SeqIO
import csv

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

def run_command(command, working_directory=None):
    print(f"Executing command: {command}")
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=working_directory
        )
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

import os

def process_sra(input_id, output_dir, ecoli_fa, threshold, include_header=True):
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

    print(f"Output files will be saved in: {working_directory}")

    # Step 1: Run fastq-dump to get reads (first 1000)
    print("\nStep 1: Running fastq-dump...")
    fastq_dump_command = f"fastq-dump --skip-technical --split-files -X 1000 {input_id} -O {working_directory}"
    if not run_command(fastq_dump_command):
        print("Error in fastq-dump command. Exiting.")
        return None

    # Debug: List files after fastq-dump
    print("\nListing files after fastq-dump:")
    run_command(f"ls -l {working_directory}")

    # Check for the presence of the files
    available_files = [f for f in os.listdir(working_directory) if f.endswith('.fastq')]
    
    if not available_files:
        print("Error: No FASTQ files found after fastq-dump.")
        return None

    # Determine which files are available
    fastq_file_1_exists = f"{input_id}_1.fastq" in available_files
    fastq_file_2_exists = f"{input_id}_2.fastq" in available_files

    # Select files to process based on availability
    if fastq_file_1_exists:
        fastq_file_1 = fastq_dump_file_1
    elif fastq_file_2_exists:
        fastq_file_1 = fastq_dump_file_2  # Fallback if _1 doesn't exist
    else:
        print("Error: Neither _1 nor _2 FASTQ files are available for processing.")
        return None

    # Proceed with processing the available file
    try:
        with open(fastq_file_1) as f:
            read_count = sum(1 for line in f) // 4
        print(f"Read count from {fastq_file_1}: {read_count}")
    except FileNotFoundError:
        print(f"Error: File {fastq_file_1} not found.")
        return None

    if read_count < 500:
        error_message = f"Error: Run {input_id} has less than 500 reads and the current sample was skipped."
        print(error_message)
        return None

    # Check if we have paired-end or single-end reads
    is_paired = os.path.exists(fastq_dump_file_1) and os.path.exists(fastq_dump_file_2)
    print(f"Is paired-end: {is_paired}")
    print(f"File 1 exists: {os.path.exists(fastq_dump_file_1)}")
    print(f"File 2 exists: {os.path.exists(fastq_dump_file_2)}")

    # Step 2: Run fastp
    print("\nStep 2: Running fastp...")
    if is_paired:
        print("Processing paired-end reads with fastp...")
        
        # Check for the existence of both files
        if os.path.exists(fastq_dump_file_1) and os.path.exists(fastq_dump_file_2):
            fastp_command = (
                f"fastp --in1 {fastq_dump_file_1} --in2 {fastq_dump_file_2} "
                f"--out1 {fastp_file_1} --out2 {fastp_file_2} --detect_adapter_for_pe"
            )
        else:
            print("One or both paired-end files not found. Exiting.")
            return None

    else:
        print("Processing single-end reads with fastp...")
        
        # Check for the first file, and fallback to the second if not found
        if os.path.exists(fastq_dump_file_1):
            fastp_command = f"fastp --in1 {fastq_dump_file_1} --out1 {fastp_file_1}"
        elif os.path.exists(fastq_dump_file_2):
            print("fastq_dump_file_1 not found, using fastq_dump_file_2...")
            fastp_command = f"fastp --in1 {fastq_dump_file_2} --out1 {fastp_file_1}"
        else:
            print("Both fastq_dump_file_1 and fastq_dump_file_2 not found. Exiting.")
            return None

    # Run the fastp command
    if not run_command(fastp_command):
        print("Error in fastp command. Exiting.")
        return None

    # Debug: List files after fastp
    print("\nListing files after fastp:")
    run_command(f"ls -l {working_directory}")

    # Step 3: Run vsearch for paired-end reads or sed for single-end reads
    if is_paired:
        print("\nStep 3: Running vsearch to merge paired-end reads...")
        vsearch_command = (
            f"vsearch --fastq_mergepairs {fastp_file_1} --reverse {fastp_file_2} --fastaout {merged_file}"
        )
        if not run_command(vsearch_command):
            print("Error in vsearch command. Exiting.")
            return None
        sed_file = merged_file
    else:
        print("\nStep 3: Running sed to convert FASTQ to FASTA...")
        sed_command = f"sed -n '1~4s/^@/>/p;2~4p' {fastp_file_1} > {sed_file}"
        if not run_command(sed_command):
            print("Error in sed command. Exiting.")
            return None

    # Debug: List files after vsearch/sed
    print("\nListing files after vsearch/sed:")
    run_command(f"ls -l {working_directory}")

    # Step 4: Align reads to reference
    print("\nStep 4: Aligning reads to reference...")
    align_command = (
        f"mafft --thread 16 --6merpair --nuc --keeplength "
        f"--addfragments {sed_file} {ecoli_fa} > {aligned_file}"
    )
    if not run_command(align_command):
        print("Error in alignment. Exiting.")
        return None

    # Step 5: Process aligned file
    print("\nStep 5: Processing aligned file...")
    hvreglocator_output = process_fasta(aligned_file, threshold, input_id, include_header=True)

    return hvreglocator_output

def run_command_simple(command, working_directory=None):
    """A simpler run_command function without debug prints."""
    try:
        subprocess.run(
            command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=working_directory
        )
        return True
    except subprocess.CalledProcessError:
        return False

def process_fasta(input_fasta, threshold, sample_id, include_header=True):
    print(f"Processing FASTA file: {input_fasta}")
    fasta_sequences = list(SeqIO.parse(input_fasta, 'fasta'))
    print(f"Total sequences read: {len(fasta_sequences)}")
    
    if len(fasta_sequences) < 2:
        print("Error: Not enough sequences in the alignment file.")
        return "Error\tInsufficient sequences for analysis\n"
    
    alignment_start = []
    alignment_end = []
    
    for seq in fasta_sequences:
        seq_str = str(seq.seq)
        # Find start position (first non-gap character)
        start_pos = next((i for i, char in enumerate(seq_str) if char != '-'), None)
        # Find end position (last non-gap character)
        end_pos = len(seq_str) - next((i for i, char in enumerate(reversed(seq_str)) if char != '-'), None) - 1
        
        if start_pos is not None and end_pos is not None:
            alignment_start.append(start_pos)
            alignment_end.append(end_pos)
    
    if not alignment_start or not alignment_end:
        print("Error: No valid alignment positions found.")
        return "Error\tNo valid alignment positions\n"
    
    # Calculate median start and end positions
    median_start = int(statistics.median(alignment_start))
    median_end = int(statistics.median(alignment_end))
    
    print(f"Median start position: {median_start}")
    print(f"Median end position: {median_end}")
    
    # Calculate coverage for all regions
    coverage_details = {}
    for region_num, region_start, region_end in regions:
        hv_length = region_end - region_start + 1
        overlap_start = max(median_start, region_start)
        overlap_end = min(median_end, region_end)
        
        if overlap_start > overlap_end:
            overlap_length = 0
        else:
            overlap_length = overlap_end - overlap_start + 1
        
        coverage = overlap_length / hv_length
        coverage_details[region_num] = coverage
    
    print("Coverage details per region:")
    for region_num in sorted(coverage_details.keys()):
        print(f"V{region_num}: {coverage_details[region_num]:.2f}")
    
    # Identify regions that meet the threshold
    covered_regions = [rn for rn, cov in coverage_details.items() if cov >= threshold]
    
    # Assign HV_region_start and HV_region_end based on covered regions
    if covered_regions:
        hv_region_start_num = min(covered_regions)
        hv_region_end_num = max(covered_regions)
        hv_region_start = f"V{hv_region_start_num}"
        hv_region_end = f"V{hv_region_end_num}"
    else:
        # Assign the region with the highest coverage if no regions meet the threshold
        max_coverage_region = max(coverage_details, key=coverage_details.get)
        max_coverage = coverage_details[max_coverage_region]
        hv_region_start = f"V{max_coverage_region}"
        hv_region_end = f"V{max_coverage_region}"
    
    # Collect warnings for regions below threshold but with some coverage
    coverage_below_threshold = [rn for rn, cov in coverage_details.items() if 0 < cov < threshold]
    
    if coverage_below_threshold:
        if len(coverage_below_threshold) == 1:
            warning_message = f"V{coverage_below_threshold[0]} below threshold of {threshold}"
        elif len(coverage_below_threshold) == 2:
            warning_message = f"V{coverage_below_threshold[0]} and V{coverage_below_threshold[1]} below threshold of {threshold}"
        else:
            regions_list = ", ".join(f"V{rn}" for rn in coverage_below_threshold[:-1])
            regions_list += f", and V{coverage_below_threshold[-1]}"
            warning_message = f"{regions_list} below threshold of {threshold}"
    else:
        warning_message = ""
    
    # Create TSV output
    headers = ["Sample_ID", "Median_Alignment_start", "Median_Alignment_end", 
               "HV_region_start", "HV_region_end", 
               "Warnings"]
    
    # Add coverage details per region as new columns with prefix Cov_
    coverage_columns = [f"Cov_V{rn}" for rn in range(1, 10)]
    headers.extend(coverage_columns)
    
    # Prepare coverage values in order from V1 to V9
    coverage_values = [f"{coverage_details.get(rn, 0):.2f}" for rn in range(1, 10)]
    
    values = [
        sample_id, 
        median_start, 
        median_end, 
        hv_region_start, 
        hv_region_end,
        warning_message
    ]
    values.extend(coverage_values)
    
    if include_header:
        results = "\t".join(headers) + "\n" + "\t".join(map(str, values)) + "\n"
    else:
        results = "\t".join(map(str, values)) + "\n"
    
    print(results)
    return results

def process_table(input_tsv, output_dir, threshold, combined_output_file):
    print(f"Processing existing TSV file: {input_tsv}")
    
    if not os.path.exists(input_tsv):
        print(f"Error: Input TSV file {input_tsv} does not exist.")
        sys.exit(1)
    
    with open(input_tsv, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames
        required_fields = ["Sample_ID", "Median_Alignment_start", "Median_Alignment_end", 
                           "HV_region_start", "HV_region_end", "Warnings"] + [f"Cov_V{rn}" for rn in range(1, 10)]
        
        if not all(field in fieldnames for field in required_fields):
            print("Error: Input TSV file is missing required columns.")
            sys.exit(1)
    
        combined_headers = ["Sample_ID", "Median_Alignment_start", "Median_Alignment_end", 
                            "HV_region_start", "HV_region_end", 
                            "Warnings"] + [f"Cov_V{rn}" for rn in range(1, 10)]
    
        with open(combined_output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=combined_headers, delimiter='\t')
            writer.writeheader()
            
            for row in reader:
                sample_id = row["Sample_ID"]
                median_start = int(row["Median_Alignment_start"])
                median_end = int(row["Median_Alignment_end"])
                
                # Extract coverage details
                coverage_details = {}
                for rn in range(1, 10):
                    cov_key = f"Cov_V{rn}"
                    try:
                        coverage_details[rn] = float(row.get(cov_key, 0))
                    except ValueError:
                        coverage_details[rn] = 0.0  # Handle non-numeric values gracefully
                
                # Identify regions that meet the new threshold
                covered_regions = [rn for rn, cov in coverage_details.items() if cov >= threshold]
                
                # Assign HV_region_start and HV_region_end based on new threshold
                if covered_regions:
                    hv_region_start_num = min(covered_regions)
                    hv_region_end_num = max(covered_regions)
                    hv_region_start = f"V{hv_region_start_num}"
                    hv_region_end = f"V{hv_region_end_num}"
                else:
                    # Assign the region with the highest coverage if no regions meet the threshold
                    max_coverage_region = max(coverage_details, key=coverage_details.get)
                    max_coverage = coverage_details[max_coverage_region]
                    hv_region_start = f"V{max_coverage_region}"
                    hv_region_end = f"V{max_coverage_region}"
                
                # Collect warnings for regions below threshold but with some coverage
                coverage_below_threshold = [rn for rn, cov in coverage_details.items() if 0 < cov < threshold]
                
                if coverage_below_threshold:
                    if len(coverage_below_threshold) == 1:
                        warning_message = f"V{coverage_below_threshold[0]} below threshold of {threshold}"
                    elif len(coverage_below_threshold) == 2:
                        warning_message = f"V{coverage_below_threshold[0]} and V{coverage_below_threshold[1]} below threshold of {threshold}"
                    else:
                        regions_list = ", ".join(f"V{rn}" for rn in coverage_below_threshold[:-1])
                        regions_list += f", and V{coverage_below_threshold[-1]}"
                        warning_message = f"{regions_list} below threshold of {threshold}"
                else:
                    warning_message = ""
                
                # Prepare the new row
                new_row = {
                    "Sample_ID": sample_id,
                    "Median_Alignment_start": median_start,
                    "Median_Alignment_end": median_end,
                    "HV_region_start": hv_region_start,
                    "HV_region_end": hv_region_end,
                    "Warnings": warning_message
                }
                
                # Add coverage details
                for rn in range(1, 10):
                    new_row[f"Cov_V{rn}"] = f"{coverage_details[rn]:.2f}"
                
                writer.writerow(new_row)
                print(f"Processed Sample_ID: {sample_id}")
    
    print(f"\nRecalculation complete. Updated results saved in: {combined_output_file}")

def process_id_list(id_list_file, output_dir, ecoli_fa, threshold):
    with open(id_list_file, 'r') as f:
        run_ids = [line.strip() for line in f if line.strip()]
    
    combined_output_file = os.path.join(output_dir, "hvreglocator_combined_output.txt")
    print(f"Combined output will be saved in: {combined_output_file}")
    
    with open(combined_output_file, 'w') as outfile:
        pass  # Just to clear the file if it exists
    
    for idx, run_id in enumerate(run_ids):
        print(f"\nProcessing run ID: {run_id}")
        hvreglocator_output = process_sra(run_id, output_dir, ecoli_fa, threshold, include_header=(idx == 0))
        if hvreglocator_output:
            with open(combined_output_file, 'a') as outfile:
                outfile.write(hvreglocator_output)
    
    print(f"\nAll results have been combined and saved in: {combined_output_file}")

def main():
    print("Starting HVRegLocator...")
    parser = argparse.ArgumentParser(
        description="HVRegLocator: Process SRA data, FASTA file, or existing TSV to locate hypervariable regions."
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Sub-command to run")
    
    # Subparser for 'sra' command
    parser_sra = subparsers.add_parser("sra", help="Process SRA Run ID(s)")
    parser_sra.add_argument(
        "-r", "--run-id", help="Single SRA Run ID for processing (required unless -l is used)"
    )
    parser_sra.add_argument(
        "-l", "--list-file", help="File containing list of SRA Run IDs to process (one per line)"
    )
    parser_sra.add_argument(
        "-e", "--ecoli", help="Path to ecoli.fa file", default=None
    )
    parser_sra.add_argument(
        "-o", "--output-dir", help="Output directory", required=True
    )
    parser_sra.add_argument(
        "-t", "--threshold", type=float, default=0.7, 
        help="Threshold for HV region coverage proportion (0-1). Default is 0.7"
    )
    
    # Subparser for 'fasta' command
    parser_fasta = subparsers.add_parser("fasta", help="Process a FASTA file")
    parser_fasta.add_argument(
        "-f", "--fasta-file", help="Input FASTA file (required)", required=True
    )
    parser_fasta.add_argument(
        "-e", "--ecoli", help="Path to ecoli.fa file", default=None
    )
    parser_fasta.add_argument(
        "-o", "--output-dir", help="Output directory", required=True
    )
    parser_fasta.add_argument(
        "-t", "--threshold", type=float, default=0.7, 
        help="Threshold for HV region coverage proportion (0-1). Default is 0.7"
    )
    
    # Subparser for 'table' command
    parser_table = subparsers.add_parser("table", help="Recalculate HV regions from an existing TSV file")
    parser_table.add_argument(
        "-i", "--input-tsv", help="Input TSV file generated by HVRegLocator (required)", required=True
    )
    parser_table.add_argument(
        "-o", "--output-dir", help="Output directory", required=True
    )
    parser_table.add_argument(
        "-t", "--threshold", type=float, default=0.7, 
        help="New threshold for HV region coverage proportion (0-1). Default is 0.7"
    )
    
    args = parser.parse_args()
    print(f"Command-line arguments: {args}")

    # Validate threshold
    if not (0.0 <= args.threshold <= 1.0):
        print("Error: Threshold must be between 0 and 1.")
        sys.exit(1)

    if args.command in ["sra", "fasta"]:
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
            process_id_list(args.list_file, args.output_dir, ecoli_fa, args.threshold)
        elif args.run_id:
            print(f"Processing single SRA run ID: {args.run_id}")
            hvreglocator_output = process_sra(args.run_id, args.output_dir, ecoli_fa, args.threshold, include_header=True)
            if hvreglocator_output:
                # Write to a single output file named after the run ID
                output_file = os.path.join(args.output_dir, f"{args.run_id}_hvreglocator_output.txt")
                with open(output_file, 'w') as f:
                    f.write(hvreglocator_output)
                print(f"Final output saved in: {output_file}")
        else:
            print("Error: The 'sra' command requires either --run-id or --list-file argument")
            parser.error("The 'sra' command requires either --run-id or --list-file argument")

    elif args.command == "fasta":
        print(f"Processing FASTA file: {args.fasta_file}")
        # Extract Sample_ID from the filename
        sample_id = os.path.splitext(os.path.basename(args.fasta_file))[0]
        hvreglocator_output = process_fasta(args.fasta_file, args.threshold, sample_id, include_header=True)
        if hvreglocator_output:
            # Write the output to a file in the output directory
            output_file = os.path.join(args.output_dir, "hvreglocator_output.txt")
            with open(output_file, 'w') as f:
                f.write(hvreglocator_output)
            print(f"Final output saved in: {output_file}")

    elif args.command == "table":
        input_tsv = os.path.abspath(args.input_tsv)
        combined_output_file = os.path.join(args.output_dir, "hvreglocator_recalculated_output.txt")
        process_table(input_tsv, args.output_dir, args.threshold, combined_output_file)
    
    print("HVRegLocator execution complete.")

if __name__ == "__main__":
    main()