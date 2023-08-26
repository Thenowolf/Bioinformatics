import itertools
import re
from datetime import datetime

blosum6 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}

blosum62 = { ("A", "A"): 5, ("T","A"): -4, ("T","T"):5, ("C","A"):-4, ("C","T"): -4, ("C","C"): 5, ("G","A"): -4, ("G","T"): -4, ("G","C"):-4, ("G","G"):5}
abeceda = {"A", "T", "C", "G"}
#abeceda = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

def makeSubstrings(seq, w):
    #dic = {i:[] for i in range(len(seq) - (w-1))}
    array = [0 for i in range(len(seq) - (w-1))] # nechci mít prázdné hodnoty
    for i in range(len(seq)):
            if(i+w <= len(seq)): # nechci mít hodnoty s délkou menší než k
                array[i] = seq[i:i+w]
                #dic[i].append(seq[i:i+w])
    return array

def compareStr(str1, str2):
    count = 0
    for i in range(len(str1)):
        try:
            count += blosum62[str2[i], str1[i]]
        except:
             #print("Nenalezeno, zkousim otocit poradi")
             count += blosum62[str1[i], str2[i]]
        #finally:
    return count

def compareStrDrop(str1, str2, dropoff, rev=False):
    if (str1 or str2) == "":
        return 0
    if rev == True:
        str1 = str1[::-1]
        str2 = str2[::-1]
    count = 0
    for i in range(len(str1)):
        if (str1[i] or str2[i]) == None:
            break
        try:
            if dropoff <= blosum62[str2[i], str1[i]]:
                count += blosum62[str2[i], str1[i]]
            else:
                break
        except:
            if dropoff <= blosum62[str1[i], str2[i]]:
             #print("Nenalezeno, zkousim otocit poradi")
                count += blosum62[str1[i], str2[i]]
            else:
                break
        #finally:
    return count

str1 = "LEH"
str2 = "QEH"

def comparePermutations(w):
    #abeceda = ["A", "B", "C"]
    arrayToCompare = list()
    for i in itertools.product(abeceda, repeat = w):
        arrayToCompare.append(''.join(i))
        #break
    return arrayToCompare

def evalseq(seq,refseq, w, T, dropoff):
    arr1 = makeSubstrings(seq, w)
    arr2 = comparePermutations(w)
    count = 0
    finalMSP = 0
    finalSeq = ""
    for i in range(len(arr1)):
        for y in range(len(arr2)):
            evaluation = compareStr(arr1[i], arr2[y])
            if evaluation >= T:
                #seqNew = seq
                #print(seq[:i])
                #print(arr1[i])
                #print(seq[i+w]:)
                #print(arr2[y])
                if refseq.find(arr2[y]) != -1:
                    print(arr2[y])
                    seqNew = seq[:i] + arr2[y] + seq[i+w:]
                    resultRight = compareStrDrop(seqNew[i:], refseq[i:], dropoff)
                #resultLeft = compareStrDrop(seqNew[:i+w], str2[:i+w], dropoff)
                    resultLeft = compareStrDrop(seqNew[:i], refseq[:i], dropoff, True)
                    MSP = resultRight + resultLeft
                    if finalMSP < MSP:
                        finalMSP = MSP
                        finalSeq = seqNew
                #seqLeft = 
                #seqRight = 
                #count += 1
                #seqNew[i:w+i] = arr1[i]
    return finalMSP, finalSeq
        #print(count)
        #count = 0

def evalseq2(seq,refseq, w, T, dropoff):
    arr1 = makeSubstrings(seq, w)
    arr2 = comparePermutations(w)
    count = 0
    finalMSP = 0
    finalSeq = ""
    finalPos = None
    for i in range(len(arr1)):
        for y in range(len(arr2)):
            evaluation = compareStr(arr1[i], arr2[y])
            if evaluation >= T:
                #seqNew = seq
                #print(seq[:i])
                #print(arr1[i])
                #print(seq[i+w]:)
                #print(arr2[y])
                if refseq.find(arr2[y]) != -1:
                    res = [k.start() for k in re.finditer(arr2[y], refseq)]
                    for z in range(len(res)):
                        print(arr2[y])
                        seqNew = seq[:i] + arr2[y] + seq[i+w:]
                        resultRight = compareStrDrop(seqNew[i:], refseq[res[z]:], dropoff)
                #resultLeft = compareStrDrop(seqNew[:i+w], str2[:i+w], dropoff)
                        resultLeft = compareStrDrop(seqNew[:i], refseq[:res[z]], dropoff, True)
                        MSP = resultRight + resultLeft
                        if finalMSP < MSP:
                            finalMSP = MSP
                            finalSeq = refseq[res[z]:res[z]+len(seq)]#seqNew
                            finalPos = res[z]
                #seqLeft = 
                #seqRight = 
                #count += 1
                #seqNew[i:w+i] = arr1[i]
    return finalMSP, finalSeq, finalPos
        #print(count)
        #count = 0

#print(blosum62["QL"])
#print(makeSubstrings("YANCLEHKMGS", 3))
#print(evalseq2("YANCLEHKMGS", "DAPCQEHKRGW", 3, 11, -2))

f = open("Bioinformatika/chr17.fa", "r")
chr17 = f.read()
chr17 = chr17.replace("\n", "")
chr17 = chr17.upper()

#print(evalseq2("GGGGACCCACACGTCT", "TTGCGCCTGCGCGGCGAGCCGGGCGCCCCCCTCCCCTCGTGGAGTCTGTGTAAAGCCGCCTGAGGCTGCGCTCTTCAGCCCCTGGGGACCCACACGTCTCGGTCGGTCTGCACCGTTTCGCCGGTCGCGACCATAATGCCCAGCCTGGCTCGGGAGGGCCGACTACCTGTTCCCCACCTCCCCCCAGCTGGCGCCCCACACTCCCAGGTGGGCTGGGCGGAGCTGTGGTTTCCGGCCCCATCTGCTCGGCCGCTGCCTTCCCGGGCTCTGGCTTAAGTTGCTTTCCCAAGAGAAAGCGAGCGGGTACCAGCGCCCCTTCCCAAGGCTTCTCCCGCCCGGGCCCAGCTGCTGTTGGTGGCGGAGGAGGGCACTGCGGGGAAGCCCCGGACGGCCCAGGTGATTCTTCTGGGCCATCACGCCCCTTCTTCGCGTGAATTCCTGCTCTTTGATGTGCACCAAAACCCTTCCCATTCCTAACGTTTCCTCCTCGTTCGCCCCGGCGTTCTTTGCACCTCCCTCGGCTCCTTTCGGGCTGTCTGTCCCTGTCTCCACTCGTCCTTCTTCGCTTCTCCCGGTTATATAACTCTTCCTCTCGCCGTGTCCTGGTTTCAACTCCACAGACTCCGCCCTGGAACAGCGCGGGGAGGGGCGGGAGAGTTGGGGACCCAGACAGATTCCGGCTGGCGGCCGGCTGGGGCACGGGGGATCCTGCAGATGGAAGACGGAGCCCCGCGGGACCCGGTGCCCCCGTCC", 11, 50, -2))
now = datetime.now()
print(evalseq2("GGGGACCCACACGTCT", chr17, 11, 50, -2))
print(datetime.now() - now)

#print(evalseq2("YANCLEHKMGS", "DAPCQEAPCGW", 3, 11, -2))

#val = 'DAPCQEHKRGW'.find('YAP')
#print(val)
#print(comparePermutations(3))
#print(compareStr(str1, str2))