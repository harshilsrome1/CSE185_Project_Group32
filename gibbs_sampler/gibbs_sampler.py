import sys
import os
import datetime
import argparse
from gibbs_algorithm import gibbs_sampler as gibbs
import time

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
    #benchmarking purposes
    start_time = time.time()
    
    parser = argparse.ArgumentParser(
        prog="gibbs_sampler",
        description="Command-line script to find shared motifs from a set of sequences"
    )

    parser.add_argument("fasta_ref", help="This is a fasta file containing all the set of sequences", metavar="FA FILE", type=str)
    parser.add_argument("kmerSize", help="Size of k-mer for motif output", type=int)
    parser.add_argument("randValue", help="Number of motifs outputted", type=int)
    parser.add_argument("repetitions", help="Number of times the procedure should repeat", type=int)

    parser.add_argument("-o", "--output", help="Output file path. Overwrites existing files by default.", metavar="FILE", type=str, default=None)
    parser.add_argument("-l", "--log", help="Log file path. Default: stdout", metavar="FILE", type=str, default=None)

    args = parser.parse_args()

    log = open(args.log, "w") if args.log else sys.stdout

    try:
        log.write("Welcome to the motif finding tool!\n")
        log.write("Start time: {}\n\n".format(datetime.datetime.now()))

        if not os.path.exists(args.fasta_ref):
            log.write("{} does not exist".format(args.fasta_ref))
            return
        
        seqs = read_fasta(args.fasta_ref)

        log.write("Using fasta: {}\n".format(args.fasta_ref))
        log.write("Starting analysis, looking for motifs...\n\n")
        
        answer = gibbs(seqs, args.kmerSize, args.randValue, args.repetitions)
        log.write("Done\n")
        log.write(str(answer))  # Ensure that 'answer' is converted to string if not already
        log.write("\n")
        
    finally:
        if log is not sys.stdout:
            log.close()  # Ensure log is always closed properly

    if args.output:
        with open(args.output, "w") as outf:
            outf.write(str(answer))  # Ensure that 'answer' is converted to string if not already

    #benchmarking results
    time_taken = time.time() - start_time
    print(time_taken)
    log.write("My program took", str(time_taken), "to run")

if __name__ == '__main__':
    main()
