"""---- module for bioinformatics----
contains functions for parsing FASTA and FASTQ files, processing strings,
and manipulaing and visualizing sequence datasets.
"""

#contains functions from Johns Hopkins University's "Algorithms for DNA Sequencing" course on Coursera

def readFASTA(filename):
    # parse and read FASTA file
    genome = ''

    with open(filename, 'r') as f:    
        for line in f:
            if line[0] != '>':
                genome += line.rstrip() 
    return genome
   



def readFASTQ(filename):
    # parse and read FASTQ file
    sequences = []
    qualities = []
    with open(filename) as f:
        while True:
            f.readline() 
            seq = f.readline().rstrip()
            f.readline()  
            qual = f.readline().rstrip() 
            
            if len(seq) == 0:  
                break
                
            sequences.append(seq) 
            qualities.append(qual)
            
    return sequences, qualities

 


def countBaseComp(seqs):  
    # returns the counts of each unique base found in a list
	count = collections.Counter()
	for seq in seqs:
    		count.update(seq)
	print(count)   
   
    
    
def phred33ToQ(qual):
    # convert ASCII code to integer base quality scores
    return ord(qual) - 33    
    


def histogram(quals):
    # create a histogram of base quality scores
    hist = [0]*len(quals[0])
    
    for qual in quals:
        for i in range(len(qual)):
            q = phred33ToQ(qual[i])
            hist[i] += q
    return hist   


    
def findGCbyPos(reads):
    # calculates GC content for each read
    gc = [0] * 100
    totals = [0] * 100
    
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
                
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
            
    return gc



def countBase(string):
    # return the number of each base present in the sequence
    counts = {'A' : 0, 'G' : 0, 'C' : 0 , 'T' : 0, 'N' : 0}
    
    for base in counts:
        counts[base] = string.count(base)
    
    return counts

    
    
def revComp(string):
    # return the reverse complement of a string
    list_str = list(string)
    list_str.reverse()
    
    baseComp = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G', 'N' : 'N'}
    list_str = [baseComp[base] for base in list_str]
      
    return ''.join(list_str)    

