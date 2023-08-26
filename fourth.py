import random
from datetime import datetime
from collections import defaultdict
from pysuffixarray.core import SuffixArray
import math


genes = ['AA', 'CC', 'GG', 'TT', 'AC', 'AG', 'AT', 'CG', 'CT', 'GT']
result3 = []
#globResultDistinct = []
f = open("Bioinformatika/chr17.fa", "r")
dna = f.read()
dna = dna.replace("\n", "")
dna = dna.replace("\r", "")

def bayesian_rule(a, b, b_a):
    a_b = (a * b)/b_a
    return a_b

def readSAM():
    f = open("Bioinformatika/bio.sam", "r")
    #bio = f.read()
    lines=f.readlines()[25:]
    result=[]
    result2=[]
    for x in lines:
        result.append(x.split('\t')[3])
        result2.append(x.split('\t')[9])
        result3.append(x.split('\t')[10])
    resultDistinct = list(set(result))
    resultDistinct.sort(reverse=True) # není dynamické, jen pro muj use case
    global globResultDistinct
    globResultDistinct = resultDistinct
    #for i in range(len(resultDistinct)):
        
    firstseq = result2[:result.count(resultDistinct[0])] # opět nedynamické
    secseq = result2[:result.count(resultDistinct[1])]

    f.close()
    return firstseq, secseq


def P_error(line, index):
    ch = result3[line][index]
    q = ord(ch)
    q = q - 33
    result = pow(10, -q/10)
    return result


def PBA(geneONE, seqchar, e):
    pba = 0
    if(geneONE == seqchar):
        pba = 1 - e # naimplementovat
    if(geneONE != seqchar):
        pba = e/3 # naimplementovat
    return pba

def PBG(gene, seqchar, e):
    pbg = 1/2 * PBA(gene[0], seqchar, e) + 1/2 * PBA(gene[1], seqchar, e)
    return pbg

def classification(refseq, gene):
    #res = "";
    if(refseq == gene[0] and refseq == gene[1]):
        return "match"
    if(refseq == gene[0] or refseq == gene[1]): # and refseq !=gene[0] or refseq != gene[1] and gene[0] == gene[1]):
        return "hetvar"
    if(gene[0] == gene[1]):
        return "homvar"
    #if(refseq == gene[0] or refseq == gene[1] and refseq !=gene[0] or refseq != gene[1] and gene[0] != gene[1]):
    return "hetvar_nonref"

def loadRefseq(seq,index):
    #dna = ""
    # f = open("Bioinformatika/chr17.fa", "r")
    # with open("Bioinformatika/chr17.fa", "r") as f:
    #     for i, line in enumerate(f):
    #         if i == 71: #int(globResultDistinct[seq]):
    #             dna = line
    #             return line
    #             break
    # dna = f.read()
    # dna = dna.replace("\n", "")
    # dna = dna.replace("\r", "")
    print(int(globResultDistinct[seq])+index)
    return dna[int(globResultDistinct[seq])+index]
    #return dna[7674791:7674800]
    #return dna

     

def PG(type):
    h = 0.001
    eerr = 0.01
    if type == "match":
       return 1 - (3*h) / 2
    if type == "hetvar":
        return h
    if type == "homvar":
        return h/2
    if type == "hetvar_nonref":
        return h * eerr



def probability(index):
    firstseq, secseq = readSAM()
    res = [0 for i in range(10)]
    #PDG = 
    for y in range(len(genes)):
        PDG = 1
        #seqchar = firstseq[i][index]
        for i in range(len(secseq)):
            e = P_error(i, index)
            seqchar = secseq[i][index] # vyřešit jinak + doplnit podmínku délky (není stejná)
            PDG *= PBG(genes[y], seqchar, e) # 4.1296e-319 a 2e-323 zde se zlomí
            #if PDG > 0:
            #print(PDG)
        #return(PDG)
        res[y] = bayesian_rule(PG(classification(loadRefseq(1, index),genes[y])), PDG, 1) # zbytečné čtu soubor xkrát
    return res
    #return PDG

def largestElement(arr):
    max = arr[0]
    y = 0
    for i in range(1,len(arr)):
        if max < arr[i]:
            max = arr[i]
            y = i
    return max, y

def findTwoBiggest(res, defLimit):
    array = res
    max1, i1 = largestElement(array)
    array.remove(max1)
    max2, i2 = largestElement(array)
    LOD = math.log(max1) - math.log(max2)
    if(LOD < defLimit):
        return genes[i1], genes[i2]
    else:
        return genes[i1]

#readSAM()
#print(loadRefseq(0,0))
#print(P_error(0, 0))
#for i in range(5,7):
    #res = all(i == 0.0 for i in probability(i))
#    print(probability(i))
#    print(i)
# for i in range(130):
#     print(probability(i))
#print(classification("A", "AC"))


#7674797 T - CC  homozigotní varianta
#7674326 C - TT homozigotní varianta
# 111-112

for i in range(111,112):
    #res = all(i == 0.0 for i in probability(i))
    print(probability(i))
    print(i)
test = [4.4926688095396735e-87, 0.0, 0.0, 4.4926688095396735e-77, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

print(findTwoBiggest(test, 1.0)) # hází error bo 0
#print(dna[7674791:7674800])
#print(dna[7674326:7674340])
print(int(globResultDistinct[1]))
