import random
from datetime import datetime
from collections import defaultdict


gen = "AGCT"


def complementary(dnasequence):
    result = ""
    for ch in dnasequence:
        if ch == "A":
            result += "T"
        if ch == "T":
            result += "A"
        if ch == "G":
            result += "C"
        if ch== "C":
            result += "G"
    return result

def findGenes():
    f = open("Bioinformatika/covid.fna", "r")
    dna = f.read()
    dna = dna.replace("\n", "")
    dna = dna.replace("\r", "")
    startgenes = 0
    genes = 0
    endgens = 0
    startcodon = "ATG"
    endcodon = {"TAA", "TAG", "TGA"}
    for i in range(0, len(dna), 3):
        if dna[i:i+3] == startcodon:
            startgenes += 1
            for i in range(i+3, len(dna), 3):
                if dna[i:i+3] in endcodon:
                    genes +=1
                    break
        #if dna[i:i+3] in endcodon:
            #endgens += 1
    print(str(startgenes) + " stargenu a  " + str(genes) + " endgenu")
#RNA vir

def translate():
    f = open("Bioinformatika/covid.fna", "r")
    dna = f.read()
    dna = dna.replace("\n", "")
    dna = dna.replace("\r", "")
    dna1 = dna[265:21553]
    table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    protein = ""
    if len(dna1) % 3 == 0:
            for i in range (0, len(dna1), 3):
                codon = dna1[i : i+3]
                protein += table[codon]
    return protein

def randSeq(length):
   chars = ["T", "A", "G", "C"]
   res = ''.join(random.choices(chars, k=length))
   print(res)
   return res

def hamming(seq1, seq2):
    count = 0
    if len(seq1) == len(seq2):
        for i in range(0, len(seq1)):
            if seq1[i] != seq2[i]:
                count += 1
    return count

def reversehamming(seq1, seq2):
    prob = []
    prob = [0 for i in range(20)] 
    count = 0
    count2 = 0
    for y in range(20):
        seq1 = randSeq(20)
        seq2 = randSeq(20)
        if len(seq1) == len(seq2):
            for i in range(0, len(seq1)):
                if seq1[i] == seq2[i]:
                    count += 1
        prob[y] = count/20
        count2 += count
        count = 0
    print(sum(prob)/len(prob))
    return count2/(20*20)

def levenshteinR(a,b):
    if len(a) == 0:
        return len(b);
    if len(b) == 0:
        return len(a)
    return min(levenshteinR(a[1:], b) + 1, levenshteinR(a, b[1:]) + 1, levenshteinR(a[1:], b[1:]) + (a[0] != b[0]))

def levenshteinD(a,b):
    a = " " + a
    b = " " + b
    t = {}
    aL = len(a)
    bL = len(b)
    for i in range(len(a)):
        t[i, 0] = i
    for j in range(bL):
        t[0, j] = j
    for j in range(1, len(b)):
        for i in range (1, len(a)):
            if a[i] == b[j]:
                t[i,j] = t[i-1, j-1]
            else:
                t[i,j] = min(t[i-1, j], t[i, j-1], t[i-1, j-1]) + 1
    return t[len(a)-1, len(b)-1]

def globalAlignment2(a,b):
    eval = defaultdict(dict)
    eval['A']['C'] = 4
    eval['A']['G'] = 2
    eval['A']['T'] = 4
    eval['A'][' '] = 8
    eval['C']['A'] = 4
    eval['C']['G'] = 4
    eval['C']['T'] = 2
    eval['C'][' '] = 8
    eval['G']['A'] = 2
    eval['G']['C'] = 4
    eval['G']['T'] = 4
    eval['G'][' '] = 8
    eval['T']['A'] = 4
    eval['T']['C'] = 2
    eval['T']['G'] = 4
    eval['T'][' '] = 8
    eval[' ']['A'] = 8
    eval[' ']['C'] = 8
    eval[' ']['G'] = 8
    eval[' ']['T'] = 8
    a = " " + a
    b = " " + b
    t = {}
    aL = len(a)
    bL = len(b)
    for i in range(len(a)):
        t[i, 0] = 8*i
    for j in range(bL):
        t[0, j] = 8*i
    for j in range(1, len(b)):
        for i in range (1, len(a)):
            if a[i] == b[j]:
                t[i,j] = t[i-1, j-1]
            else:
                t[i,j] = min(t[i-1, j] + eval[a[i-1]][' '], t[i, j-1] + eval[' '][b[j-1]], t[i-1, j-1] + eval[a[i-1]][b[j-1]]) #+ eval[a[i]][b[j]]
    return t[len(a)-1, len(b)-1]



def globalAlignment(x, y):
    alphabet = ["A", "C", "G", "T"] 
    score = [[0, 4, 2, 4, 8],
        [4, 0, 4, 2, 8],
        [2, 4, 0, 4, 8],
        [4, 2, 4, 0, 8],
        [8, 8, 8, 8, 0]]
    D = []
    traceback = []
    for i in range(len(x)+1):
        D.append([0]* (len(y)+1))
    for i in range(len(x)+1):
        traceback.append([0]* (len(y)+1))

    for i in range(1, len(x)+1):
        #test = score[2][3]
        #hodnota = score[alphabet.index(x[i-1])][-1]
        #D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1]
        D[i][0] = 8*i
        traceback[i][0] = "up"
    for i in range(len(y)+1):
        #hodnota = score[-1][alphabet.index(y[i-1])]
        #D[0][i] = D[0][i-1]+ score[-1][alphabet.index(y[i-1])]
        D[0][i] = 8*i
        traceback[0][i] = "left"

    traceback[0][0] = "done"

    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1]+ score[-1][alphabet.index(y[j-1])]
            distVer = D[i-1][j]+ score[alphabet.index(x[i-1])][-1]
            if x[i-1] == y[j-1]:
                #D[i,j] = D[i-1][j-1]
                #traceback[i-1][j-i] = "diag"
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            D[i][j] = min(distHor, distVer, distDiag)
            if D[i][j] == distHor:
                traceback[i][j] = "left"
            if D[i][j] == distVer:
                traceback[i][j] = "up"
            if D[i][j] == distDiag:
                traceback[i][j] = "diag"

    i = len(x)
    j = len(y)
    cigar = ""
    while(i > 0 or j>0):
        if traceback[i][j] == "diag":
            cigar += "M"
            i = i - 1
            j = j - 1
        elif traceback[i][j] == "left":
            cigar += "I"
            j = j-1
        elif traceback[i][j] == "up":
            cigar += "I"
            i = i-1
        elif traceback[i][j] == "done":
            break
    return D[-1][-1], traceback,cigar

def localAlignment(x, y):
    alphabet = ["A", "C", "G", "T"] 
    score = [[2, -4, -4, -4, -6],
            [-4, 2, -4, -4, -6],
            [-4, -4, 2, -4, -6],
            [-4, -4, -4, 2, -6],
            [-6, -6, -6, -6, 0]]
    D = []
    traceback = []
    for i in range(len(x)+1):
        D.append([0]* (len(y)+1))
    for i in range(len(x)+1):
        traceback.append([0]* (len(y)+1))

    for i in range(1, len(x)+1):
        #test = score[2][3]
        #hodnota = score[alphabet.index(x[i-1])][-1]
        #D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1]
        D[i][0] = 0*i
        traceback[i][0] = "up"
    for i in range(len(y)+1):
        #hodnota = score[-1][alphabet.index(y[i-1])]
        #D[0][i] = D[0][i-1]+ score[-1][alphabet.index(y[i-1])]
        D[0][i] = 0*i
        traceback[0][i] = "left"

    traceback[0][0] = "done"

    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1]+ score[-1][alphabet.index(y[j-1])]
            distVer = D[i-1][j]+ score[alphabet.index(x[i-1])][-1]
            if x[i-1] == y[j-1]:
                #D[i,j] = D[i-1][j-1]
                #traceback[i-1][j-i] = "diag"
                distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]
            else:
                distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            D[i][j] = max(distHor, distVer, distDiag, 0)
            if i == 1 and j == 1:
                tmp = D[i][j]
                fp = i
                sp = j
            if tmp < D[i-1][j-1]:
                tmp = D[i-1][j-1] # drobny zadrhle xd
                fp = i-1
                sp = j-1
            if D[i][j] == distHor:
                traceback[i][j] = "left"
            if D[i][j] == distVer:
                traceback[i][j] = "up"
            if D[i][j] == distDiag:
                traceback[i][j] = "diag"

    i = fp
    j = sp
    cigar = ""
    while(i > 0 or j>0):
        if traceback[i][j] == "diag":
            cigar += "M"
            i = i - 1
            j = j - 1
        elif traceback[i][j] == "left":
            cigar += "I"
            j = j-1
        elif traceback[i][j] == "up":
            cigar += "I"
            i = i-1
        elif D[i][j] == 0:
            break
    return D[-1][-1], traceback,cigar

#for i in range(10,20,1):
#print(reversehamming(randSeq(20), randSeq(20)))

#print(levenshteinR(randSeq(8), randSeq(8)))
#print(levenshteinR("ACTCGGCT", "ACATATAC"))
#print(lev("ACTCGGCT", "ACATATAC"))
#print(levenshteinD("ACTCGGCT", "ACATATAC"))

#a = randSeq(1100)
#b = randSeq(1100)

#now = datetime.now()
#print(levenshteinR(a, b)) # 10
#print(datetime.now() - now)
#now = datetime.now()
#print(levenshteinD(a, b)) #1100 cca
#print(datetime.now() - now)


eval = defaultdict(dict)
eval['A']['C'] = 4
eval['A']['G'] = 2
eval['A']['T'] = 4
eval['A'][' '] = 8
eval['C']['A'] = 4
eval['C']['G'] = 4
eval['C']['T'] = 2
eval['C'][' '] = 8
eval['G']['A'] = 2
eval['G']['C'] = 4
eval['G']['T'] = 4
eval['G'][' '] = 8
eval['T']['A'] = 4
eval['T']['C'] = 2
eval['T']['G'] = 4
eval['T'][' '] = 4
eval[' ']['A'] = 8
eval[' ']['C'] = 8
eval[' ']['G'] = 8
eval[' ']['T'] = 8

#x = 'A'
#y = 'C'
#print(eval[x][y])
#https://stackoverflow.com/questions/53784533/global-alignment-sequence-function
print(globalAlignment("TACGTCAGC", "TATGTCATGC"))
#print(localAlignment("GGTATGCTGGCGCTA", "TATATGCGGCGTTT"))



#print(hamming(randSeq(7), randSeq(7)))
#print(gen)
#print(complementary(gen))
#findGenes()
#print(randSeq(10))
#print(translate())