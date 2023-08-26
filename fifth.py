import itertools



def longestSubstringFinder(string1, string2, l):
    answer = ""
    len1, len2 = len(string1), len(string2)
    for i in range(len1):
        match = ""
        for j in range(len2):
            if (i + j < len1 and string1[i + j] == string2[j]):
                match += string2[j]
            else:
                if (len(match) > len(answer)): answer = match
                match = ""
    if len(answer) >= l:
        return answer
    else: return ""

def overlap(a, b, l):
    start = 0
    indexes = []
    lena = len(a)
    i = 0
    while True:
        if i == 0:
            start = a.find(b[:1])
            if start == -1:
                return 0
        if b.startswith(a[start:]):
            if(l <= len(a) - start):
                return len(a) - start
            else:
                return 0
        start = start + 1
        #lena = lena -1
        #i = i + 1
        #i = i - 1
        #indexes[i] = start
        i= i+1

def overlap2(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], 10)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def pickmaxoverlap(reads, l):
    maxoverlap = 0
    amax = None
    bmax = None
    #for i in range(len(reads)-1):
    for a,b in itertools.permutations(reads, 2):
        lenoverlap = overlap(a, b, l)
        if lenoverlap > maxoverlap:
            maxoverlap = lenoverlap
            amax = a
            bmax = b
    return maxoverlap, amax, bmax

def shortest_super_string(reads, l):
    while True:
        maxoverlap, a, b = pickmaxoverlap(reads, l)
        if maxoverlap > 0 :
            reads.remove(a)
            reads.remove(b)
            reads.append(a + b[maxoverlap:])
        else:
            break
    return ''.join(reads)








def shortestString(array, l): # fix chybu pro l větší než k-mer
    length = len(array)
    merge = array[0]
    for i in range(1,(len(array))):
        common = longestSubstringFinder(merge, array[i], l)
        if common == "":
            merge += array[i]
        else:
            merge += "".join(array[i].split(common)) # magie
    return merge

def makeKmers(seq, k):
    array = [0 for i in range(len(seq) - (k-1))] # nechci mít prázdné hodnoty
    for i in range(len(seq)):
        try:
            if(i+k <= len(seq)): # nechci mít hodnoty s délkou menší než k
                array[i] = seq[i:i+k]
        except:
            print("Mimo index, nevadííí")
    return array


array= ["AAA", "AAB", "ABB", "BBB", "BBA"]
array2 = ["long_lon", "ng_long_", "_long_lo", "g_long_t", "ong_long", "g_long_l", "ong_time", "a_long_l", "_long_ti", "long_tim"]
#print(shortestString(array,2))
#print(overlap(array[0], array[4], 1))
#print(shortest_super_string(array2, 3))
#print(makeKmers("CTTTCCACTTGATAAGAGGTCCCAAGACTTAGTACCTGGAGGGTGAAATATTCTCCATCCAGTGGTTTCTTCTTTGGCTGGGGAGAGGAGCTGGTGTTGTTGGGCAGTGCTAGGAAAGAGGCAAGGAAAGGTGATAAAAGTGAATCTGAGG", 4))

print(shortest_super_string(makeKmers("CTTTCCACTTGATAFAGAGGTCCCAAGACTTAGTACCTGGAGGGTGAAATATTCTCCATCCAGTGGTTTCTTCTTTGGCTGGGGAGAGGAGCTGGTGTTGTTGGGCAGTGCTAGGAAAGAGGCAAGGAAAGGTGATAAAAGTGAATCTGAGG", 8), 1))
