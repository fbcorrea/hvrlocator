import joblib
import numpy as np
import argparse
import statistics
import subprocess
import os
import sys
from Bio import SeqIO
import csv
import time
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from statistics import mean, median, stdev
from scipy.stats import skew, kurtosis


# Define hypervariable regions globally
regions = [
    (1, 0, 96),
    (2, 97, 299),
    (3, 300, 479),
    (4, 480, 811),
    (5, 812, 885),
    (6, 886, 1065),
    (7, 1066, 1180),
    (8, 1181, 1372),
    (9, 1373, 1600)
]

FEATURE_COLUMNS = [
    '1_5_count','1_5_mean','1_5_median','1_5_std','1_5_min','1_5_max','1_5_25th','1_5_75th',
    '6_10_count','6_10_mean','6_10_median','6_10_std','6_10_min','6_10_max','6_10_25th','6_10_75th'
]

def get_stats(vals, prefix):
    if len(vals) == 0:
        return {f"{prefix}_{key}": np.nan for key in [
            "count", "mean", "median", "std", "min", "max", "25th", "75th"
        ]}
    return {
        f"{prefix}_count": len(vals),
        f"{prefix}_mean": np.mean(vals),
        f"{prefix}_median": np.median(vals),
        f"{prefix}_std": np.std(vals),
        f"{prefix}_min": np.min(vals),
        f"{prefix}_max": np.max(vals),
        f"{prefix}_25th": np.percentile(vals, 25),
        f"{prefix}_75th": np.percentile(vals, 75)
    }



def run_command(command, working_directory=None):
    print(f"Executing command: {command}")
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            cwd=working_directory
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


def run_command_with_retry(command, working_directory=None, max_retries=3):
    for attempt in range(max_retries):
        print(f"Attempt {attempt + 1} of {max_retries}: Executing command: {command}")
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                cwd=working_directory
            )
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                print("Command executed successfully.")
                return True, None
            else:
                print(f"Error executing command. Return code: {process.returncode}")
                print(f"Error message: {stderr.decode('utf-8')}")
        except Exception as e:
            print(f"Exception occurred: {str(e)}")

        if attempt < max_retries - 1:
            print("Retrying in 3 seconds...")
            time.sleep(3)

    return False, stderr.decode('utf-8') if 'stderr' in locals() else str(e)

def process_sra(input_id, output_dir, ecoli_fa, threshold, model_path=None, include_header=True, report_dir=None):
    print(f"Processing SRA data for ID: {input_id}")
    working_directory = os.path.join(output_dir, input_id)
    os.makedirs(working_directory, exist_ok=True)

    fastq1 = os.path.join(working_directory, f"{input_id}_1.fastq")
    fastq2 = os.path.join(working_directory, f"{input_id}_2.fastq")
    fastp1 = os.path.join(working_directory, f"{input_id}_fastp_1.fastq")
    fastp2 = os.path.join(working_directory, f"{input_id}_fastp_2.fastq")
    merged = os.path.join(working_directory, f"{input_id}_merged.fasta")
    processed = os.path.join(working_directory, f"{input_id}_processed.fasta")
    aligned = os.path.join(working_directory, f"{input_id}_aligned.fasta")

    if report_dir:
        os.makedirs(report_dir, exist_ok=True)
        fastp_html = os.path.join(report_dir, f"{input_id}_fastp.html")
        fastp_json = os.path.join(report_dir, f"{input_id}_fastp.json")
    else:
        fastp_html = os.path.join(working_directory, f"{input_id}_fastp.html")
        fastp_json = os.path.join(working_directory, f"{input_id}_fastp.json")

    cmd = f"fastq-dump --skip-technical --split-files -X 1000 {input_id} -O {working_directory}"
    success, err = run_command_with_retry(cmd, working_directory)
    if not success:
        return False, f"Error in fastq-dump: {err}"

    files = [f for f in os.listdir(working_directory) if f.endswith('.fastq')]
    if not files:
        return False, f"No FASTQ files found for {input_id}"

    paired = os.path.exists(fastq1) and os.path.exists(fastq2)

    if paired:
        primary_fastq = fastq1
    else:
        if os.path.exists(fastq1):
            primary_fastq = fastq1
        elif os.path.exists(fastq2):
            primary_fastq = fastq2
        else:
            return False, f"FASTQ files not found for {input_id}"

    with open(primary_fastq) as fh:
        read_count = sum(1 for _ in fh) // 4
    if read_count < 500:
        return False, f"Run {input_id} has less than 500 reads. Skipped."

    if paired:
        cmd = (
            f"fastp --in1 {fastq1} --in2 {fastq2} "
            f"--out1 {fastp1} --out2 {fastp2} "
            f"--detect_adapter_for_pe --html {fastp_html} --json {fastp_json}"
        )
    else:
        cmd = (
            f"fastp --in1 {primary_fastq} --out1 {fastp1} "
            f"--html {fastp_html} --json {fastp_json}"
        )
    if not run_command(cmd):
        return False, f"Error in fastp for {input_id}"

    primer_flag = 'NA'
    primer_score = 'NA'
    if model_path and os.path.exists(model_path):
        try:
            rf_model = joblib.load(model_path)
            s1, s2 = [], []
            with open(fastp1, 'rt') as fh:
                for rec in SeqIO.parse(fh, 'fastq'):
                    q = rec.letter_annotations['phred_quality']
                    if len(q) >= 10:
                        s1.extend(q[:5])
                        s2.extend(q[5:10])
            a1, a2 = np.array(s1), np.array(s2)
            m1 = get_stats(a1, prefix='1_5')
            m2 = get_stats(a2, prefix='6_10')
            feat = [m1[c] for c in FEATURE_COLUMNS[:8]] + [m2[c] for c in FEATURE_COLUMNS[8:]]
            proba = rf_model.predict_proba([feat])[0][1]
            primer_flag = 'TRUE' if proba > 0.5 else 'FALSE'
            primer_score = f"{proba:.3f}"
        except Exception as e:
            print(f"Model error: {e}")

    if paired:
        cmd = f"vsearch --fastq_mergepairs {fastp1} --reverse {fastp2} --fastaout {merged}"
        if not run_command(cmd):
            return False, f"Error in vsearch for {input_id}"
        seq_input = merged
    else:
        cmd = f"sed -n '1~4s/^@/>/p;2~4p' {fastp1} > {processed}"
        if not run_command(cmd):
            return False, f"Error in sed for {input_id}"
        seq_input = processed

    cmd = (
        f"mafft --thread 16 --6merpair --nuc --keeplength "
        f"--addfragments {seq_input} {ecoli_fa} > {aligned}"
    )
    if not run_command(cmd):
        return False, f"Error in alignment for {input_id}"

    output_line = process_fasta(aligned, threshold, input_id, include_header, primer_flag, primer_score)
    return True, output_line

def process_fasta(input_fasta, threshold, sample_id, include_header=True,model_path=None, primer_flag='NA', primer_score='NA'):
     if model_path and os.path.exists(model_path) and (
        input_fasta.endswith(".fastq") or input_fasta.endswith(".fq") or input_fasta.endswith(".fastq.gz")
    ):
        try:
            rf_model = joblib.load(model_path)
            s1, s2 = [], []
            # Use FastqGeneralIterator for speed; handles plain FASTQ (not gz)
            # If gzipped, you'd need to open with gzip.open
            opener = open
            if input_fasta.endswith(".gz"):
                import gzip
                opener = gzip.open

            with opener(input_fasta, "rt") as fh:
                for header, seq, qual in FastqGeneralIterator(fh):
                    q = [ord(c) - 33 for c in qual]
                    if len(q) >= 10:
                        s1.extend(q[:5])
                        s2.extend(q[5:10])

            if s1 and s2:
                a1, a2 = np.array(s1), np.array(s2)
                m1 = get_stats(a1, prefix='1_5')
                m2 = get_stats(a2, prefix='6_10')
                feat = [m1[c] for c in FEATURE_COLUMNS[:8]] + [m2[c] for c in FEATURE_COLUMNS[8:]]
                proba = rf_model.predict_proba([feat])[0][1]
                primer_flag = 'TRUE' if proba > 0.5 else 'FALSE'
                primer_score = f"{proba:.3f}"
            else:
                print("Warning: Not enough quality scores to compute primer features; leaving primer fields as NA.")
        except Exception as e:
            print(f"Model error: {e}")
    
    seqs = list(SeqIO.parse(input_fasta, 'fasta'))
    if len(seqs) < 2:
        return "Error\tInsufficient sequences\n"

    starts, ends = [], []
    for rec in seqs:
        s = str(rec.seq)
        start = next((i for i,c in enumerate(s) if c != '-'), None)
        end = len(s) - next((i for i,c in enumerate(reversed(s)) if c != '-'), None) - 1
        if start is not None and end is not None:
            starts.append(start)
            ends.append(end)

    med_s = int(statistics.median(starts))
    med_e = int(statistics.median(ends))
    avg_s = round(statistics.mean(starts))
    avg_e = round(statistics.mean(ends))
    min_s = min(starts)
    max_e = max(ends)

    covs = {}
    for num, r0, r1 in regions:
        length = r1 - r0 + 1
        overlap = max(0, min(med_e, r1) - max(med_s, r0) + 1)
        covs[num] = overlap / length

    covered = [n for n,c in covs.items() if c >= threshold]
    if covered:
        hv_s = f"V{min(covered)}"; hv_e = f"V{max(covered)}"
    else:
        m = max(covs, key=covs.get)
        hv_s = hv_e = f"V{m}"

    below = [n for n,c in covs.items() if 0 < c < threshold]
    if below:
        if len(below) == 1:
            warn = f"V{below[0]} below threshold of {threshold}"
        elif len(below) == 2:
            warn = f"V{below[0]} and V{below[1]} below threshold of {threshold}"
        else:
            lst = ", ".join(f"V{n}" for n in below[:-1])
            warn = f"{lst}, and V{below[-1]} below threshold of {threshold}"
    else:
        warn = ""

    ps_num = next((n for n,r0,r1 in regions if r0 <= med_s <= r1), None)
    pe_num = next((n for n,r0,r1 in regions if r0 <= med_e <= r1), None)
    ps = f"V{ps_num}" if ps_num else ""; pe = f"V{pe_num}" if pe_num else ""
    cov_ps = covs.get(ps_num, 0.0); cov_pe = covs.get(pe_num, 0.0)

    headers = [
       "Sample_ID", "Primer_Presence", "Score_Primer_Presence",
        "Min_Alignment_start", "Max_Alignment_end",
        "Average_Alignment_start", "Average_Alignment_end",
        "Median_Alignment_start", "Median_Alignment_end",
        "Predicted_HV_region_start", "Predicted_HV_region_end",
        "Coverage_based_HV_region_start", "Coverage_based_HV_region_end",
        "Coverage_HV_region_start", "Coverage_HV_region_end",
        "Warnings"
    ] + [f"Cov_V{n}" for n in range(1, 10)] 
    values = [
        sample_id, primer_flag, primer_score,
        min_s, max_e,avg_s, avg_e, med_s, med_e,
        ps, pe, hv_s, hv_e,
        f"{cov_ps:.2f}", f"{cov_pe:.2f}", warn
    ] + [f"{covs[n]:.2f}" for n in range(1, 10)]

    if include_header:
        return "\t".join(headers) + "\n" + "\t".join(map(str, values)) + "\n"
    else:
        return "\t".join(map(str, values)) + "\n"


def process_table(input_tsv, output_dir, threshold, combined_output_file):
    print(f"Processing existing TSV file: {input_tsv}")
    if not os.path.exists(input_tsv):
        print(f"Error: Input TSV {input_tsv} not found.")
        sys.exit(1)
    with open(input_tsv, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        # already includes predicted and coverage columns
        headers = reader.fieldnames if reader.fieldnames else []
        with open(combined_output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=headers, delimiter='\t')
            writer.writeheader()
            for row in reader:
                writer.writerow(row)
    print(f"Recalculation complete. Saved: {combined_output_file}")

def process_id_list(id_list_file, output_dir, ecoli_fa, threshold, model_path=None, report_dir=None):
    runs = [r.strip() for r in open(id_list_file) if r.strip()]
    combined = os.path.join(output_dir, 'hvreglocator_combined_output.txt')
    errors = os.path.join(output_dir, 'hvreglocator_error_report.txt')
    os.makedirs(output_dir, exist_ok=True)
    header = [
        "Sample_ID", "Primer_Presence", "Score_Primer_Presence",
        "Min_Alignment_start", "Max_Alignment_end",
	"Average_Alignment_start", "Average_Alignment_end",
        "Median_Alignment_start", "Median_Alignment_end",
        "Predicted_HV_region_start", "Predicted_HV_region_end",
	"Coverage_based_HV_region_start", "Coverage_based_HV_region_end", 
        "Coverage_HV_region_start", "Coverage_HV_region_end",
        "Warnings"
    ] + [f"Cov_V{n}" for n in range(1, 10)]
    with open(combined, 'w') as co, open(errors, 'w') as er:
        co.write("\t".join(header) + "\n")
        er.write("Sample_ID\tError\n")
    retries = []
    for rid in runs:
        ok, res = process_sra(rid, output_dir, ecoli_fa, threshold, model_path=model_path, include_header=False, report_dir=report_dir)
        if ok:
            with open(combined, 'a') as co:
                co.write(res)
        else:
            retries.append((rid, res))
    for rid, res in retries:
        ok, res2 = process_sra(rid, output_dir, ecoli_fa, threshold, model_path=model_path, include_header=False, report_dir=report_dir)
        if ok:
            with open(combined, 'a') as co:
                co.write(res2)
        else:
            with open(errors, 'a') as er:
                er.write(f"{rid}\t{res2}\n")
    print(f"All combined: {combined}")
    print(f"Errors: {errors}")

def main():
    parser = argparse.ArgumentParser(description="HVRegLocator: locate hypervariable regions")
    sub = parser.add_subparsers(dest='command', required=True)

    sra_p = sub.add_parser('sra')
    sra_p.add_argument('-r', '--run-id')
    sra_p.add_argument('-l', '--list-file')
    sra_p.add_argument('-e', '--ecoli', default=None)
    sra_p.add_argument('-o', '--output-dir', required=True)
    sra_p.add_argument('-t', '--threshold', type=float, default=0.6)
    sra_p.add_argument('-m', '--model', help='Path to trained RF model (.pkl)')
    sra_p.add_argument('--report-dir')

    fasta_p = sub.add_parser('fasta')
    fasta_p.add_argument('-f', '--fasta-file', required=True)
    fasta_p.add_argument('-e', '--ecoli', default=None)
    fasta_p.add_argument('-o', '--output-dir', required=True)
    fasta_p.add_argument('-t', '--threshold', type=float, default=0.6)
    fasta_p.add_argument('-m', '--model', help='Path to trained RF model (.pkl)')

    table_p = sub.add_parser('table')
    table_p.add_argument('-i', '--input-tsv', required=True)
    table_p.add_argument('-o', '--output-dir', required=True)
    table_p.add_argument('-t', '--threshold', type=float, default=0.6)

    args = parser.parse_args()

    if args.command in ['sra', 'fasta']:
        ecoli_fa = args.ecoli or os.path.join(os.path.dirname(__file__), 'ecoli.fa')
        if not os.path.exists(ecoli_fa):
            print(f"E. coli reference not found: {ecoli_fa}")
            sys.exit(1)
        os.makedirs(args.output_dir, exist_ok=True)

    if args.command == 'sra':
        if args.list_file:
            process_id_list(args.list_file, args.output_dir, ecoli_fa, args.threshold, model_path=args.model, report_dir=args.report_dir)
        elif args.run_id:
            ok, res = process_sra(args.run_id, args.output_dir, ecoli_fa, args.threshold, model_path=args.model, include_header=True, report_dir=args.report_dir)
            if ok:
                outf = os.path.join(args.output_dir, f"{args.run_id}_hvreglocator_output.txt")
                with open(outf, 'w') as f:
                    f.write(res)
                print(f"Output saved: {outf}")
            else:
                print(res)
        else:
            parser.error('sra requires --run-id or --list-file')

    elif args.command == 'fasta':
        sid = os.path.splitext(os.path.basename(args.fasta_file))[0]
        res = process_fasta(args.fasta_file, args.threshold, sid, include_header=True, model_path=args.model)
        outf = os.path.join(args.output_dir, 'hvreglocator_output.txt')
        with open(outf, 'w') as f:
            f.write(res)
        print(f"Output saved: {outf}")

    elif args.command == 'table':
        inp = os.path.abspath(args.input_tsv)
        outp = os.path.join(args.output_dir, 'hvreglocator_recalculated_output.txt')
        process_table(inp, args.output_dir, args.threshold, outp)

if __name__ == '__main__':
    main()
