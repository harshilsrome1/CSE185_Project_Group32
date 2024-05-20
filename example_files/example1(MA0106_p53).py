#this file contains an example using a p53 gene in homo sapiens 
#we pulled the 17 binding sites of the transcription factor from the JASPAR database, and are running gibbs sampler on it to return the top motifs, based on some parameters we have put

dna4_test = []
#this is the set of 17 sequences from JASPAR
dna1 = "gcggggacattagggaCGGGACATGCCCGGGCATGTct"
dna2 = "ctaaggcccccatcCGGGACATGCCCGGGCATGTt"
dna3 = "gaAGGGACATGCCCGGGCATGTcggaaacgata"
dna4 = "tgaatattgaggcttcACGGACATGTCCGGGCATGTct"
dna5 = "acgaaatccttcATGGACATGCCCGGGCATGTcg"
dna6 = "gcccactaccaatgagCGGGACATGTCCGGACATGTct"
dna7 = "gacccacacctaatgtCAAGACATGCCCGGGCATGTct"
dna8 = "ACAAACATGCCCGGGCATGTtcctgctgtggacct"
dna9 = "ctgcCCAAACATGTCCGGGCATGTccgtgacgtgg"
dna10 = "gaaaaatagatggtgGGGAACATGCCCGGGCATGTct"
dna11 = "cgcaaacgaTAGAACATGCCCGGGCATGTccaggg"
dna12 = "cggtactctgcatgaTAGGACATGTCCGGACATGTtc"
dna13 = "agcgaactgtaagtgcGGGGGCATGCCCGGGCATGTct"
dna14 = "ACAAACATGCCCGGGCATGCcccatgtga"
dna15 = "ggacatgcCCGGGCATGTCTGTAGTTCCctctac"
dna16 = "gggcatgcCCGGGCATGTTCTTGGAATAacaacct"
dna17 = "ggacatgcCCGGGCATGCCCAAAATAAGaccatta"

#just to run gibbs sampler, we used capital letters for all nucleotides, so I use the upper function to ensure all strings are all uppercase
dna1 = dna1.upper()
dna2 = dna2.upper()
dna3 = dna3.upper()
dna4 = dna4.upper()
dna5 = dna5.upper()
dna6 = dna6.upper()
dna7 = dna7.upper()
dna8 = dna8.upper()
dna9 = dna9.upper()
dna10 = dna10.upper()
dna11 = dna11.upper()
dna12 = dna12.upper()
dna13 = dna13.upper()
dna14 = dna14.upper()
dna15 = dna15.upper()
dna16 = dna16.upper()
dna17 = dna17.upper()

#we add all the final sequences into one list to run gibbs sampler on
dna4_test.append(dna1)
dna4_test.append(dna2)
dna4_test.append(dna3)
dna4_test.append(dna4)
dna4_test.append(dna5)
dna4_test.append(dna6)
dna4_test.append(dna7)
dna4_test.append(dna8)
dna4_test.append(dna9)
dna4_test.append(dna10)
dna4_test.append(dna11)
dna4_test.append(dna12)
dna4_test.append(dna13)
dna4_test.append(dna14)
dna4_test.append(dna15)
dna4_test.append(dna16)
dna4_test.append(dna17)

#this is where we run the gibbs sampler
#we used some test parameters, let me explain what each is 
#the first parameter is the k-value, in this case we are looking for 6-mer motifs
#the next parameter is how many motifs we want to output, in this case we set it to 8 so it will return the top 8 motifs
#the final parameter is how many times the procedure should repeat in looking for the set of best k-mer motifs, in this case we set it to 100
#  for the final parameter, the tradeoff is that a large number will be more guarenteed to be correct, but it will take much longer to run 
answer4 = gibbs_sampler(dna4_test,6,8,100)
print(answer4)
