import numpy as np
from Bio import SeqIO
import array as arr


reads = np.empty(0,dtype=str)
qualities = np.empty(0, dtype=int)
n = len(list(SeqIO.parse("test.txt", "fastq")))
qualities = np.empty(n, dtype = list)
index = 0

for record in SeqIO.parse("test.txt", "fastq"):

    reads = np.append(reads, str(record.seq))
    quality = np.array(record.letter_annotations["phred_quality"])
    # qualities = np.append(qualities, quality)
    qualities[index] = quality
    index += 1

print(reads[0])
print(qualities)
