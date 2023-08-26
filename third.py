import random
from datetime import datetime
from collections import defaultdict
from pysuffixarray.core import SuffixArray
import sais


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


# Quick sort in Python

# function to find the partition position
def partition(array, low, high):

  # choose the rightmost element as pivot
  temp = array[high]
  array[high] = array[low]
  array[low] = temp

  pivot = array[high]

  # pointer for greater element
  i = low

  # traverse through all elements
  # compare each element with pivot
  for j in range(low, high):
    if array[j] <= pivot:
      # if element smaller than pivot is found
      # swap it with the greater element pointed by i
      #i = i + 1

      # swapping element at i with element at j
      (array[i], array[j]) = (array[j], array[i])
      i = i + 1

  # swap the pivot element with the greater element specified by i
  (array[i], array[high]) = (array[high], array[i])

  # return the position from where partition is done
  return i

# function to perform quicksort
def quickSort(array, low, high):
  if low < high:

    # find pivot element such that
    # element smaller than pivot are on the left
    # element greater than pivot are on the right
    pi = partition(array, low, high)

    # recursive call on the left of pivot
    quickSort(array, low, pi - 1)

    # recursive call on the right of pivot
    quickSort(array, pi + 1, high)

def alphabetSort(arr):
    for i in range(len(arr)):
        for j in range(i+1, len(arr)):
            #if compareTo(arr[i], arr[j]) > 0:
            #if arr[i].compareTo(arr[j]) > 0 :
            #if ord(arr[i]) - ord(arr[j]) > 0:
            if arr[i] > arr[j]:
                temp = arr[i]
                arr[i] = arr[j]
                arr[j] = temp            
    return arr

def count_sort_letters(array, size, col, base, max_len):
  """ Helper routine for performing a count sort based upon column col """
  output   = [0] * size
  count    = [0] * (base + 1) # One addition cell to account for dummy letter
  min_base = ord('a') - 1 # subtract one too allow for dummy character

  for item in array: # generate Counts
    # get column letter if within string, else use dummy position of 0
    letter = ord(item[col]) - min_base if col < len(item) else 0
    count[letter] += 1

  for i in range(len(count)-1):   # Accumulate counts
      count[i + 1] += count[i]

  for item in reversed(array):
    # Get index of current letter of item at index col in count array
    letter = ord(item[col]) - min_base if col < len(item) else 0
    output[count[letter] - 1] = item
    count[letter] -= 1

  return output

def radix_sort_letters(array, max_col = None):
  """ Main sorting routine """
  if not max_col:
    max_col = len(max(array, key = len)) # edit to max length

  for col in range(max_col-1, -1, -1): # max_len-1, max_len-2, ...0
    array = count_sort_letters(array, len(array), col, 4, max_col)

  return array
 
def suffixArray(s):
    sa = [0 for i in range(len(s))]
    arr = [0 for i in range(len(s))] 
    for i in range(len(s)):
        arr[i] = s[i:]
        #suffix= sorted([s[i:]])
    #arr = alphabetSort(arr)
    #arr.sort()
    quickSort(arr, 0, len(arr)-1)
    #arr = radix_sort_letters(arr)
    for i in range(len(s)):
        sa[i] = len(s)-len(arr[i])
    return sa

def binarySearchSA(t, sa, p):
    if len(t) == 1: return 1 #aby to zbytečně nepočítalo
    l = 0 
    r = len(sa) 
    while True:
        c = (l + r) // 2 # prostředek + zaokrouhlení
        plt = True # assume p < T[sa[c]:] until proven otherwise
        i = 0
        while i < len(p) and sa[c]+i < len(t): # dokud neporovnávám prázdné znaky
            if p[i] < t[sa[c]+i]:
                break
            elif p[i] > t[sa[c]+i]:
                plt = False
                break
            i += 1 # tied so far
        if plt:
            if c == (l + 1):
                return c # není už kde jinde půlit => našel jsem
            r = c
        else:
            if c == (r - 1):
                return r # r protože jsem ho už jednou našel, ale nevěděl jsem o tom
            l = c

def bisect_left(a, x, text):
    pocetit = 0
    l = 0
    r = len(a)
    #Binary search
    while l < r:
        mid = (l+r)//2
        #porovnávání stringů - jednoduše, dokud nenalezne první nesrovnalost
        if text[a[mid]:] == x:
            return a[mid]
        elif text[a[mid]:] < x: 
            pocetit += 1 
            l = mid+1 # aby nenastala situace, že budu porovnávat stejný prvek -> spatny vysledek
            #l = mid # alternativa - chtěl jsem vyzkoušet
            #if l == r - 1:
                #l = l + 1
        else: 
        	r = mid
    if text[a[l]:].startswith(x): 
        # i suppose text[a[lo]:a[lo]+len(x)] == x could be a faster check
        #print(pocetit)
        return a[l]
    #lo has index of first match
    else:
        return "not found"

def BWT(seq, sufarr):
    L = ""
    F = ""
    for i in range(len(seq)):
        F += seq[sufarr[i]]
        if sufarr[i] == 0:
            L += '$'
        if sufarr[i] > 0:
            L += seq[sufarr[i]-1]
    #print(L)
    #print(F)
    return L

def rank(L):
    unique_chars = list(set(L))
    len_unique_chars = len(unique_chars)
    arr = [0 for i in range(len_unique_chars)]
    for i in range(len_unique_chars):
        arr[i] = [0 for i in range(len(L))]
    for i in range(len(L)):
        for y in range(len_unique_chars):
            if L[i] == unique_chars[y] :
                arr[y][i] = arr[y][i-1] + 1
            elif L[i] != unique_chars[y]:
                arr[y][i] += arr[y][i-1]
    return arr

def rankBwt(bw):
    tots	= dict()
    ranks	=	[]
    for	c	in	bw:
        if	c	not in	tots:
            tots[c]	= 0
        ranks.append(tots[c])
        tots[c]	+= 1
    return	 tots

def firstCol(tots):
	first	=	{}
	totc	= 0
	for	c, count in sorted(tots.items()):
		first[c]	=	(totc,	totc	+	count)
		totc	+=	count
	return	first

def rank(bwt, c, k): # je potřeba dopočítat pozici, abych se pak mohl podívat, kde se nachází v F
        return bwt[:k].count(c)

def rank_lt(F, c): # nalezne začátek charu
    if c in F:
        return F[c][0]
    else:
        return None

def search(sa, bwt, pat):
        F = firstCol(rankBwt(bwt))
        L = bwt
        begin = 0
        end = len(L)
        for c in pat[::-1]:
            offset = rank_lt(F, c)
            if offset is None:
                begin, end = None, None
                break
            begin = offset + rank(bwt, c, begin) # at vím kde se dívat v F == SA
            end   = offset + rank(bwt, c, end)
            #begin = rank_lt(F, c, 0)
            #end = rank_lt(F, c, 1)
            if begin >= end: # no results
                begin, end = None, None
                break
        # print('[bwt] (begin, end)', begin, end)
        match = []
        if begin is not None and end is not None:
            for i in range(begin, end):
                match.append((sa[i], sa[i] + len(pat))) # vrátí výskyt všech podstringu
        return match

def readSAM():
    f = open("Bioinformatika/bio.sam", "r")
    #bio = f.read()
    lines=f.readlines()[25:]
    result=[]
    result2=[]
    for x in lines:
        result.append(x.split('\t')[9])
        result2.append(x.split('\t')[10])
    f.close()
    return result,result2
    #bio = chr17.upper()

        




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
#print(globalAlignment("TACGTCAGC", "TATGTCATGC"))
#print(localAlignment("GGTATGCTGGCGCTA", "TATATGCGGCGTTT"))
#print(suffixArray("abaaba$"))
#print(binarySearchSA("ab$", suffixArray("ab$"), "$"))
#print(bisect_left(suffixArray("abaaba$"),"baa","abaaba$"))

#f = open("Bioinformatika/chr17.fa", "r")
#chr17 = f.read()
#chr17 = chr17.upper()
#now = datetime.now()
#seq = sais.Sequence(bytes(chr17, 'utf-8'))
#print(search(seq.suffix_array[:len(seq)], BWT(chr17,seq.suffix_array[:len(seq)]) , "CAGGACTGCTCGAGC"))
#seq.suffix_array[:len(seq)] # 0:00:13.859984 bez printu
#print(datetime.now() - now)

print(readSAM())

#seq = sais.Sequence(bytes(chr17, 'utf-8'))
#seq = sais.Sequence(bytes("abaaba$", 'utf-8'))
#print(seq)
#now = datetime.now()
#print(bisect_left(seq.suffix_array[:len(seq)],"CAGGACTGCTCGAGC",chr17))
#print(bisect_left(seq.suffix_array[:len(seq)],"abaa","abaaba$"))
#print(datetime.now() - now)
#print(rank(BWT("abaaba$", suffixArray("abaaba$"))))
#print(rankBwt(BWT("abaaba$", suffixArray("abaaba$"))))
#print(firstCol(rankBwt(BWT("abaaba$", suffixArray("abaaba$")))))
#print(search(suffixArray("abaaba$"), BWT("abaaba$",suffixArray("abaaba$")) , "aba"))
#print(rankBwt("abaaba$"))
#print(chr17[7822115:7823115])
#print(suffixArray(chr17))
#now = datetime.now()
#print(bisect_left(SuffixArray(chr17),"CTTTCCACTTGATAAGAGGTCCCAAGACTTAGTACCTGGAGGGTGAAATATTCTCCATCCAGTGGTTTCTTCTTTGGCTGGGGAGAGGAGCTGGTGTTGTTGGGCAGTGCTAGGAAAGAGGCAAGGAAAGGTGATAAAAGTGAATCTGAGG",chr17))
#print(datetime.now() - now)
#print(hamming(randSeq(7), randSeq(7)))
#print(gen)
#print(complementary(gen))
#findGenes()
#print(randSeq(10))
#print(translate())