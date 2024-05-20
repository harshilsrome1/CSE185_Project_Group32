import sys 
import random
import numpy as np
import math 
import argparse
import datetime
import pandas as pd

def main():

    # Command Line Stuff
    parser = argparse.ArgumentParser(
        prog = "gibbs_sampler",
        description = "Command-line script to find shared motifs from a set of sequences"
    )

    # input arguments 
    parser.add_argument("fasta_ref",\
                        help="This is a fasta file containing all the set of sequences", \
                        metavar="FA FILE", type = str, required = True) 
    
    parser.add_argument("kmerSize", help = "size of k-mer we want in motifs output", type = int, required = True)
    parser.add_argument("randValue", help = "number of motifs outputted", type = int, required = True)
    parser.add_argument("repetitions", help = "amount of times we want procedure to repeat", type = int, required = True)

    # Output
    parser.add_argument("-o", "--output", help="Output file path."\
                        "This file will be overwritten if it already exists."\
                        "Default: stdout", metavar="FILE", type = str, required=False)
    parser.add_argument("-l", "--log", help = "write log to file. Default: stdout", metavar = "FILE", type = str, required = False)


    args = parser.parse_args()

    # output file set-up
    if args.output is None:
        outf = sys.stdout
    else:
        outf = open(args.output, "w")

    # set up log file
    if args.log is None:
        log = sys.stdout
    else:
        log = open(args.log, "w")

    log.write("Welcome to the motif finding tool!\n")
    log.write("Start time: ")
    log.write(str(datetime.datetime.now()))
    log.write("\n\n")

    # leading input files

    log.write("Loading input files:\n")

    # load fasta file
    if args.fasta_ref is not None:
        if not os.path.exists(args.fasta_ref):
            myutils.ERROR("{fasta} does not exist".format(fasta = args.fasta_ref))
        try:
            reffasta = pyfaidx.Fasta(args.fasta_ref)
        except Exception:
            myutils.ERROR("please check fasta file format - see README.md for details")
        log.write("Using fasta: {fasta}".format(fasta = args.fasta_ref))
        log.write("\n")
    else:
        myutils.ERROR("please specify a fasta file")
    
def gibbs_sampler(dna: list[str], k: int, t: int, n: int) -> list[str]:
    bestAnswers = []
    for eachdna in dna:
        number = rand(len(eachdna)-k)
        first = eachdna[number:number+k]
        bestAnswers.append(first)
    for j in range(1000):
        motifs = []
        bestMotifs = []
        for eachdna in dna:
            number = rand(len(eachdna)-k)
            first = eachdna[number:number+k]
            motifs.append(first)
        bestMotifs = motifs
        matrix = []
        for a in range(n+1):
            probDis = []
            modMotifs = []
            numberi = rand(t)
            rowControl = 0
            for a in range(len(motifs)-1):
                if a == numberi:
                    rowControl += 1
                modMotifs.append(motifs[rowControl])
                rowControl += 1
            matrix = profile(modMotifs)
            text = dna[numberi]
            indices = []
            for cc in range(len(text)-k):
                pattern = text[cc:cc+k]
                indices.append(cc)
                prob = stringProb(pattern,matrix)
                probDis.append(prob)
            randIndex = random.choices(indices,probDis)
            randomIndex = randIndex[0]
            kmer = text[randomIndex:randomIndex+k]
            motifs[numberi] = kmer
            if(score(motifs) < score(bestMotifs)):
                bestMotifs = motifs
        if(score(bestMotifs) < score(bestAnswers)):
            bestAnswers = bestMotifs;
    return(bestAnswers)

def profile(motifs) -> list[dict[str, float]]:
    answer = []
    size = len(motifs)
    eachSize = len(motifs[0])
    matrix = []
    nucleotides = ["A", "C", "G", "T"]
    sized = size+4
    for z in range(eachSize):
        eachDict = {}
        for y in nucleotides:
            eachDict[y] = 1/sized
        matrix.append(eachDict)
    for d in range(size):
        for q in range(eachSize):
            if(motifs[d][q] == "A"):
                matrix[q]["A"] += 1/sized
            if(motifs[d][q] == "C"):
                matrix[q]["C"] += 1/sized
            if(motifs[d][q] == "G"):
                matrix[q]["G"] += 1/sized
            if(motifs[d][q] == "T"):
                matrix[q]["T"] += 1/sized
    return matrix
    
def rand(n: int) -> int: 
    return(random.randint(0, n-1))

def score(motifs):
    score = 0
    size = len(motifs[0])
    num = len(motifs)
    
    for m in range(size):
        a = 0
        c = 0
        g = 0
        t = 0
        for x in range(len(motifs)):
            if motifs[x][m] == "A":
                a += 1
            if motifs[x][m] == "C":
                c += 1
            if motifs[x][m] == "G":
                g += 1
            if motifs[x][m] == "T":
                t += 1
        score += (len(motifs)-max(a,c,g,t))
    return score

def profile_most_probable_kmer(text: str, k: int, profile: list[dict[str, float]]) -> str:
    n = len(text)
    percent = 1
    maxpercent = 0
    for i in range(n-k+1):
        pattern = text[i:i+k]
        percent = stringProb(pattern,profile)
        if percent > maxpercent:
            maxpercent = percent
            answer = pattern
    return answer

def stringProb(pattern: str, profile: list[dict[str, float]]):
    percent = 1
    value = 0
    for l in range(len(pattern)):
        if pattern[l] == "A":
            value = profile[l].get("A")
            percent *= value
        if pattern[l] == "C":
            value = profile[l].get("C")
            percent *= value
        if pattern[l] == "G":
            value = profile[l].get("G")
            percent *= value
        if pattern[l] == "T":
            value = profile[l].get("T")
            percent *= value
    return percent

#test1
a = "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"
b = "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"
c = "TAGTACCGAGACCGAAAGAAGTATACAGGCGT"
d = "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"
e = "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
dna = []
dna.append(a)
dna.append(b)
dna.append(c)
dna.append(d)
dna.append(e)
answer = gibbs_sampler(dna,8,5,100)
print(answer) 
# Output for answer = ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']

#test2
dna2 = []
dna2.append("AATTGGCACATCATTATCGATAACGATTCGCCGCATTGCC")
dna2.append("GGTTAACATCGAATAACTGACACCTGCTCTGGCACCGCTC")
dna2.append("AATTGGCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT")
dna2.append("GGTTAATGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG")
dna2.append("AATTGGACGGCAACTACGGTTACAACGCAGCAAGAATATT")
dna2.append("GGTTAACTGTTGTTGCTAACACCGTTAAGCGACGGCAACT")
dna2.append("AATTGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG")
dna2.append("GGTTAAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA")
answer2 = gibbs_sampler(dna2,6,8,100)
print(answer2)
# Output for answer2 = ['AATTGG', 'AACTGA', 'AATTGG', 'AAGTGT', 'AATTGG', 'AACTGT', 'AATTGG', 'AAAAGG']

#test3
dna3 = []
dna3.append("GCACATCATTAAACGATTCGCCGCATTGCCTCGATTAACC")
dna3.append("TCATAACTGACACCTGCTCTGGCACCGCTCATCCAAGGCC")
dna3.append("AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTAACC")
dna3.append("AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAAGGCC")
dna3.append("AACCGGACGGCAACTACGGTTACAACGCAGCAAGTTAACC")
dna3.append("AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGAAGGCC")
dna3.append("AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTTTAACC")
dna3.append("AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCC")
answer3 = gibbs_sampler(dna3,6,8,100)
print(answer3)
# Output for answer3 = ['ATTAAC', 'CATAAC', 'CTTAAC', 'GTTATG', 'GTTAAC', 'GCTAAC', 'TTTAAC', 'GCTTTC']
