import sys
from Bio import SeqIO

# Usage "align_pos.py FASTA"

# Input argument
input_fasta = sys.argv[1]

sequences = []
for record in SeqIO.parse(input_fasta, "fasta"):
    sequences.append(record)
sequences = sequences[1:] # Dropping first sequence (ref full-length ecoli 16s)

sequence_counts = {}
# Count the occurrences of each sequence (ignoring /1 or /2)
for record in sequences:
    sequence_id = record.id.rsplit('/', 1)[0]  # Remove the /1 or /2 part
    sequence_counts[sequence_id] = sequence_counts.get(sequence_id, 0) + 1
ids_to_keep = []

for k,v in sequence_counts.items():
    if v == 2:
        ids_to_keep.append(k)

# Filter the sequences
filtered_sequences = []
for record in sequences:
    if record.id.rsplit('/', 1)[0] in ids_to_keep:
        filtered_sequences.append(record)


# Define variables
counts = [] # store counts for each sequence
countsmat = [] # store all counts for all sequences in matrix
count = 1

# Read each fasta sequence and count the number of spacers ("-")
for record in sequences:
	# Count the number of spacers at the beginning of the sequence
	iteration = 0
	for a, b in zip(list(record.seq), list(record.seq)[1:]):
		if(a=="-"):
			if a==b:
				count += 1
			else:
				counts.append((count))
				count = 1
		else:
			if iteration == 0:
				counts.append(1)
#			counts.append(1)
			break
		iteration =+ 1
	# Count the number of spacers at the end of the sequence (reversed)
	iteration = 0
	for a, b in zip(list(reversed(record.seq)), list(reversed(record.seq))[1:]):
		if(a=="-"):
			if a==b:
				count += 1
			else:
				counts.append((len(record.seq)-count))
				count = 1
		else:
			if iteration == 0:
				counts.append(1)
#			counts.append(1)
			break
		iteration += 1
	#print(counts)
	countsmat.append(counts)
	count=1
	counts=[]


# Calculate average of all start / end positions at alignment
column_start = []; # isolate column 0
column_end = []; # isolate column 1

for row in countsmat:
  column_start.append(row[0])

#start_pos_avg = int(round(sum(column_start)/len(column_start)))
for row in countsmat:
  column_end.append(row[1])
#end_pos_avg = int(round(sum(column_end)/len(column_end)))

import numpy as np
from scipy import stats
np.set_printoptions(threshold=sys.maxsize)

column_start = np.array(column_start)
column_end = np.array(column_end)
#print(np.sort(column_start))
#print(np.sort(column_end))

mode_start = stats.mode(column_start)
mode_end = stats.mode(column_end)

start_pos_avg = mode_start.mode.astype(int)
end_pos_avg = mode_end.mode.astype(int)

# Figure out spanning region (for Paired-end)
region_start = None
region_end = None

# Forward primer
if start_pos_avg <= 69:
	region_start = 1
elif start_pos_avg >= 99 and start_pos_avg <= 157:
	region_start = 2
elif start_pos_avg >= 227 and start_pos_avg <= 440:
	region_start = 3
elif start_pos_avg >= 500 and start_pos_avg <= 590:
	region_start = 4
elif start_pos_avg >= 650 and start_pos_avg <= 828:
	region_start = 5
elif start_pos_avg >= 857 and start_pos_avg <= 1000:
	region_start = 6
elif start_pos_avg >= 1036 and start_pos_avg <= 1119:
	region_start = 7
elif start_pos_avg >= 1158 and start_pos_avg <= 1243:
	region_start = 8
elif start_pos_avg >= 1295 and start_pos_avg <= 1435:
	region_start = 9

# Reverse primer
if end_pos_avg >= 1465:
	region_end = 9
elif end_pos_avg <= 1435 and end_pos_avg >= 1295:
	region_end = 8
elif end_pos_avg <= 1243 and end_pos_avg >= 1158:
	region_end = 7
elif end_pos_avg <= 1119 and end_pos_avg >= 1036:
	region_end = 6
elif end_pos_avg <= 1000 and end_pos_avg >= 857:
	region_end = 5
elif end_pos_avg <= 828 and end_pos_avg >= 650:
	region_end = 4
elif end_pos_avg <= 590 and end_pos_avg >= 500:
	region_end = 3
elif end_pos_avg <= 440 and end_pos_avg >= 227:
	region_end = 2
elif end_pos_avg <= 157 and end_pos_avg >= 99:
	region_end = 1


print("Alignment start:",start_pos_avg)
print("Alignment end:", end_pos_avg)
print("Hypervariable region span: ","V"+str(region_start)+"-V"+str(region_end))
