"""
Authors: Robert Schmidt, Kelby Kies, Shatabdi Sen, Parnal Joshi
3/9/2020
HW 3
"""


import numpy as np
from Bio import SeqIO
import math

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
    return eta_out



def GenS(reads):
    """
    Generate First guess at true sequences
    :param reads:
    :return:
    """
    s_out = np.empty(0,dtype="str")
    for i in range(0, K):
        x=np.random.randint(0, len(reads))
        s_out=np.append(s_out, reads[x])
    s_out = np.ndarray.tolist(s_out)
    return s_out


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
    e_ik = math.log(pie_k)
    for p in range(0,len(read)):
        r_ip = read[p]
        s_kp = s_k[p]
        a = code[s_kp]
        b = code[r_ip]
        q_ip = q_i[p]

        if r_ip == s_kp:
            e_ik += math.log((1-10**(-q_ip/10)))

        else:
            e_ik += math.log((10**(-q_ip/10)*eta[a][b]))

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


    E_t = np.ndarray((len(reads), K))
    for i in range(0, len(reads)):
        esum=0
        E_i = np.zeros(K)
        read = reads[i]

        quality = qualities[i]


        for k in range(0, K):

            # Numerator
            E_i[k] = EstepNumerator(pie[k], eta, str(s[k]),  str(read), list(quality))

        # print("Num", E_i)
        max_value = np.max(E_i)

        E_i = np.exp(E_i)
        E_i = E_i - max_value

        esum = np.sum(E_i)
        E_it = E_i/esum
        E_t[i] = E_it

    return np.asarray(E_t)


def UpdatePi(E):
    """
    M step update function for pi
    :param E:e_ik(t) matrix of expectations of P(Z=k|Data and Parameters)
    :return:vector of new pis
    """
    outPie = np.zeros(K)
    for k in range(0, K):
        relevantEs = E[:, k]
        relevantEsum = np.sum(relevantEs)
        outPie[k] = relevantEsum/len(reads)
    return outPie


def UpdateEta(reads, s_k, E):
    """
    :param reads:
    :param s_k:
    :param E:
    :return:
    """

    eta_out = np.zeros((4,4))
    for i in range(0,len(reads)):
        eta_num = np.zeros((4, 4))

        relevantE = E[i][0]
        read = reads[i]

        for p in range(0, len(read)):
            readBase = code[read[p]]
            trueBase = code[s_k[p]]
            if readBase == trueBase:
                continue
            else:
                eta_num[trueBase][readBase] += 1

        eta_out += eta_num

    sum_of_etas = eta_out.sum(axis=1)[:,None]
    for i in range(0, 4):
        if sum_of_etas[i] == 0:
            continue
        else:
            eta_out[i] = eta_out[i]/sum_of_etas[i]

    return eta_out


def UpdateS(reads, qualities, E, eta):
    sNew = np.empty(0,dtype="str")
    l = len(reads[1])
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
        sNew = np.append(sNew, s_kNew)
    return sNew


def Convergence(pieNew, pieOld, etaNew, etaOld):
    """

    :param pieNew:
    :param pieOld:
    :param etaNew:
    :param etaOld:
    :return:
    """
    global pie_is_con
    global eta_is_con
    # Test Pie
    if max(pieNew - pieOld) < 10**(-6):
        print("")
        print("Max Pie Difference",max(pieNew - pieOld))
        print("")
        pie_is_con = True

    if np.max(etaNew - etaOld) < 10**(-6):
        print("")
        print("Eta Difference", np.max(etaNew - etaOld))
        print("")
        eta_is_con = True

    return pie_is_con and eta_is_con




## Main ##

# Intiailize Best guess at iteration 0
K=12
code={"A":0, "C":1, "G":2, "T":3 }
# Import Fastq
reads = []
qualities = []
pie = np.random.dirichlet(np.ones(K))

for record in SeqIO.parse("real_amplicons.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])
eta = GenEta()
s = GenS(reads)

# Run EM algorithm
pie_is_con = False
eta_is_con = False

pieOld = [0]*K
pieNew = pie
etaOld = [[0,0,0,0]]*4
etaNew = eta
iteration = 0
while not Convergence(pieNew, pieOld, etaNew, etaOld):


    print("E1")
    E = Estep(pieNew, etaNew, reads, qualities, s)

    # M1Step:
    pieTemp = UpdatePi(E)
    s = UpdateS(reads, qualities, E, etaNew)

    # Estep
    E = Estep(pieTemp, etaNew, reads, qualities, s)

    # M2Step:
    pieOld = pieNew
    pieNew = UpdatePi(E)
    etaOld = etaNew

    etaNew = UpdateEta(reads, s[0], E)
    iteration += 1


# Output True sequences and proportions in FASTA
fh = open('Output_FASTA_real.fasta', 'w')
for i in range(0, K):
    fh.write(">TrueSeq:")
    fh.write(str(i+1))
    fh.write(" ")
    fh.write(str(pie[i]))
    fh.write("\n")
    fh.write(s[i])
    fh.write("\n")

fh.close()
