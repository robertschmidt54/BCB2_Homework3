import numpy as np
from Bio import SeqIO

# K=2
# code={"A":0, "C":1, "G":2, "T":3 }

# reads = np.empty(0,dtype='str')
#
# qualities=np.empty(0, dtype=int)
# for record in SeqIO.parse("TestData.fastq", "fastq"):
#     print(record.seq)
#     reads = np.append(reads, str(record.seq))
#     temp=np.array(record.letter_annotations["phred_quality"])
#     qualities=np.append(qualities, temp, 1)


def GenEta():
    """
    Generate first guess at Eta Matrix
    :return:
    """
    eta_out = np.zeros((4,4))
    for i in range(0,4):
        eta_temp = np.random.dirichlet(np.ones(3), size=1)
        if i == 3:
            eta_out[i]=np.append(eta_temp, 0)
        else:
            eta_out[i]=np.insert(eta_temp, i, 0, 1)
    eta_out = np.ndarray.tolist(eta_out)
    return eta_out
# eta=GenEta()


def GenS(reads):
    """
    Generate First guess at true sequences
    :param reads:
    :return:
    """
    s_out=np.empty(0,dtype="str")
    for i in range(0, K):
        x=np.random.randint(0, len(reads))
        s_out=np.append(s_out, reads[x])
    s_out = np.ndarray.tolist(s_out)
    return s_out

# s=GenS(reads)


#Test Data:
# K=2
# reads=["ATTA", "ATGC", "GGGG", "GGCG", "ATAT", "TGGT", "AACT", "ATTT", "ATGC"]
# qualities=[[40, 20, 10,1], [32, 22, 11,1], [42, 42, 42, 42],[42, 42, 42, 2],[42, 42, 42, 3],[42, 42, 42, 4],[42, 42, 42, 5],[42, 42, 42, 6], [42, 42, 42,7]]
# pie=[0.3, 0.3, 0.4]
# eta=[[0, 0.3, 0.3, 0.4],
#      [0.3, 0, 0.3, 0.4],
#      [0.3, 0.3, 0, 0.4],
#      [0.3, 0.3, 0.4, 0]]
# s=["TTTT", "CCCC"]


def EstepNumerator(pie_k, eta, s_k, read, q_i):
    """
    :param pie: vector of pis
    :param eta: matrix of etas
    :param sk: true sequence string
    :param reads: vector of i reads string
    :param quailty: vector of quailty scores for reads numeric
    :param n: number of reads
    :param K: number of clusters/distributions/species/what ever.
    :return: a single number e_ik numerator.
    """
    e_ik=1
    for p in range(0,len(read)):
        r_ip=read[p]
        s_kp=s_k[p]
        a=code[s_kp]
        b=code[r_ip]
        q_ip=q_i[p]

        if r_ip == s_kp:
            e_ik*=(1-10**(-q_ip/10))*pie_k

        else:
            e_ik*=(10**(-q_ip/10)*eta[a][b])*pie_k

    return e_ik

def Estep(pie, eta, reads, qualities, s):
    """

    :param pie: vector of pis
    :param eta: Matrix of eta values
    :param reads: vector of reads
    :param qualities: vector of quality scores
    :param s: vector of true sequences
    :return: matrix containing all e_ik
    """
    # print("E STEP:")
    E_t = np.ndarray((len(reads), K))
    for i in range(0, len(reads)):
        esum=0
        logsum=0
        E_i = np.zeros(K)
        # E_i = np.empty(0)
        read=reads[i]
        quality=qualities[i]
        for k in range(0,K):

            E_i[k] = EstepNumerator(pie[k], eta, str(s[k]),  str(read), list(quality))
            # print("E Step: k", k)
            # print(E_i[k])
            esum += E_i[k]

        E_it=E_i/esum
        logsum += np.log(esum)

        E_t[i] = E_it

    return np.asarray(E_t)


def UpdatePi(E):
    """
    M step update function for pi
    :param E:e_ik(t) matrix of expectations of P(Z=k|Data and Parameters)
    :return:vector of new pis
    """
    outPie=np.zeros(K)
    for k in range(0, K):
        relevantEs=E[:,k]
        relevantEsum=np.sum(relevantEs)
        outPie[k]=relevantEsum/len(reads)
    return outPie


def UpdateEta(reads, s_k, E):
    """

    :param reads:
    :param s_k:
    :param E:
    :return:
    """

    eta_out=np.zeros((4,4))
    # print("Reads: Update Eta", reads)
    # print("S: Update Eta",s_k)
    # print("E: Update Eta",E)
    # print(len(reads))
    for i in range(0,len(reads)):
        eta_num = np.zeros((4, 4))
        #relevantE=E[i]
        relevantE=E[i]
        print("relevant E", relevantE)
        read=reads[i]
        for p in range(0,len(read)):
            readBase=code[read[p]]
            trueBase=code[s_k[p]]
            if readBase == trueBase:
                continue
            else:
                eta_num[trueBase][readBase] += 1

        eta_out+=(eta_num*relevantE)
        # print("Numerator", eta_num)
        # print("relevantE", relevantE)
    # print(eta_out)
    # print(eta_out[:])
    # print(np.sum(eta_out, axis=1))

    eta_out=eta_out[:]/np.sum(eta_out,axis=1)

    return eta_out


def UpdateS(reads, qualities, E, eta):
    sNew = np.empty(0,dtype="str")
    l=len(reads[1])
    for k in range(0, K):
        s_kNew = ""
        for p in range(0, l):
            sumLikihood = np.zeros(4)
            for nuc in ["A", "C", "G", "T"]:

                for i in range(0, len(reads)):
                    read = reads[i]
                    quality = qualities[i]
                    readBase = code[read[p]]
                    if readBase == code[nuc]:
                        sumLikihood[code[nuc]] += np.log((1 - 10 ** (-quality[p] / 10)))*E[i][k]
                    else:
                        sumLikihood[code[nuc]] += (np.log(10 ** (-quality[p] / 10)) + np.log(eta[code[nuc]][readBase]))*E[i][k]
            s_kp = list(code.keys())[list(code.values()).index(np.argmax(sumLikihood))]
            s_kNew += s_kp
        sNew=np.append(sNew, s_kNew)
    return sNew

## Main ##

# Intiailize Best guess at iteration 0
K=2
code={"A":0, "C":1, "G":2, "T":3 }
# Import Fastq
reads = []
qualities = []
pie = np.ndarray.tolist(np.random.dirichlet(np.ones(K)))
eta=GenEta()
s=GenS(reads)
for record in SeqIO.parse("test.txt", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])

# Run EM algorithm
E="none"

for t in range(0, 10):
    # print("Iteration "+str(t)+":\n")
    # print("E = ", E)
    # print("pie = ",pie)
    # print("eta = ",eta)
    # print("s= ", s)
    # Estep:
    E=Estep(pie, eta, reads, qualities, s)
    # print("First E Step:", E)
    # M1Step:
    pie=UpdatePi(E)
    s=UpdateS(reads, qualities, E, eta)
    # Estep
    E=Estep(pie, eta, reads, qualities, s)
    # print("Second E Step:", E)
    # M2Step:
    pie=UpdatePi(E)
    # print(E)
    eta=UpdateEta(reads, s[1], E)

# Output True sequences and proportions in FASTA
fh = open('Output_FASTA.fasta', 'w')
for i in range(0,K):
    fh.write(">TrueSeq:")
    fh.write(str(i))
    fh.write(" ")
    fh.write(str(pie[i]))
    fh.write("\n")
    fh.write(s[i])
    fh.write("\n")

fh.close()