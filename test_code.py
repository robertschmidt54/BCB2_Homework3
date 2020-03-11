import numpy as np
from Bio import SeqIO
# import array as arr
#
#
# # reads = np.empty(0,dtype=str)
# # qualities = np.empty(0, dtype=int)
# n = len(list(SeqIO.parse("test.txt", "fastq")))
# # qualities = np.empty(n, dtype = list)
# reads = []
# qualities = []
# for record in SeqIO.parse("test.txt", "fastq"):
#
#     # reads = np.append(reads, str(record.seq))
#     # quality = np.array(record.letter_annotations["phred_quality"])
#     # qualities = np.append(qualities, quality)
#     # qualities[index] = quality
#
#     reads.append(str(record.seq))
#     # qualities.append(rec)
#     qualities.append(record.letter_annotations["phred_quality"])
#     # print(type(record.letter_annotations["phred_quality"]))
#
#
# print(type(reads))
# reads2=["ATTA", "ATGC", "GGGG", "GGCG", "ATAT", "TGGT", "AACT", "ATTT", "ATGC"]
# print(type(reads2))
# # print(type(qualities))
#
# s=["TTTT", "CCCC"]
# # print(type(s))
#
# qualities=[[40, 20, 10,1], [32, 22, 11,1], [42, 42, 42, 42],[42, 42, 42, 2],[42, 42, 42, 3],[42, 42, 42, 4],[42, 42, 42, 5],[42, 42, 42, 6], [42, 42, 42,7]]
# # print(type(qualities))
#
# pie1=[0.3, 0.3, 0.4]
# # print(type(pie1[0]))
# K = 2
# pie2 =  np.ndarray.tolist(np.random.dirichlet(np.ones(K),size=1))
# # print(type(pie2[0][0]))
# def GenEta():
#     eta_out = np.zeros((4,4))
#     for i in range(0,4):
#         eta_temp = np.random.dirichlet(np.ones(3), size=1)
#         if i == 3:
#             eta_out[i]=np.append(eta_temp, 0)
#         else:
#             eta_out[i]=np.insert(eta_temp, i, 0, 1)
#     eta_out = np.ndarray.tolist(eta_out)
#     return eta_out
# eta=GenEta()
# # print(type(eta))
#
# # print(type(eta))
# def GenS(reads):
#     s_out=np.empty(0,dtype="str")
#     for i in range(0, K):
#         x=np.random.randint(0, len(reads))
#         s_out=np.append(s_out, reads[x])
#     s_out = np.ndarray.tolist(s_out)
#     return s_out
#
# s=GenS(reads)
# # print("S", type(s))
#
#
# # Output FASTA
# fh  = open('Output_FASTA.fasta', 'w')
# for i in range(0,K):
#     fh.write(">TrueSeq:")
#     fh.write(str(i))
#     fh.write(" ")
#     fh.write(str(pie1[i]))
#     fh.write("\n")
#     fh.write(s[i])
#     fh.write("\n")
#
# fh.close()




e = np.random.rand(2,2)
print(e)
print(e/e.sum(axis=1)[:,None])