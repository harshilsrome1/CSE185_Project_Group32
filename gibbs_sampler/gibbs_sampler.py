import sys 
import random
import numpy as np
import math 
import argparse
import datetime
import pandas as pd
from gibbs_algorithm import gibbs_sampler as gibbs
import utils as utils

def read_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        seqs, current_seq = [], []
        for line in file:
            if line.startswith('>'):
                if current_seq:
                    seqs.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            seqs.append(''.join(current_seq))
    return seqs
    
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
    
    parser.add_argument("kmerSize", help = "size of k-mer we want in motifs output", type = int)
    parser.add_argument("randValue", help = "number of motifs outputted", type = int)
    parser.add_argument("repetitions", help = "amount of times we want procedure to repeat", type = int)

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
            utils.ERROR("{fasta} does not exist".format(fasta = args.fasta_ref))
        try:
            reffasta = read_fasta(args.fasta_ref)
        except Exception:
            utils.ERROR("please check fasta file format - see README.md for details")
        log.write("Using fasta: {fasta}".format(fasta = args.fasta_ref))
        log.write("\n")
    else:
        utils.ERROR("please specify a fasta file")

    if args.kmerSize is not None:
        kSize = kmerSize
    else:
        utils.ERROR("kmer size is not specified")

    if args.randValue is not None:
        numMotifs = randValue
    else:
        utils.ERROR("number of motifs is not specified")

    if args.repetitions is not None:
        amtRuns = repetitions
    else:
        utils.ERROR("number of runs of algorithm is not specified")

    # run the algorithm 
    log.write("Starting motif enrichment analysis...\n\n")

    answer = gibbs.gibbs_sampler(reffasta,kmerSize,numMotifs,amtRuns)
    log.write("Done\n\n")

    log.write(answer)

if __name__ == '__main__':
    main()

#test1
"""
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
"""
