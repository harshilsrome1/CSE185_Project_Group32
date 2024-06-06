import random 
from tqdm import tqdm

#this is the main gibbs sampler algorithm, it takes in 4 inputs and returns a list with the set of motifs we are looking for
def gibbs_sampler(dna: list[str], k: int, t: int, n: int) -> list[str]:
    bestAnswers = [] #this stores the best set of motifs to return 
    #this for loop creates the first set of random motifs 
    for eachdna in dna:
        number = rand(len(eachdna)-k)
        first = eachdna[number:number+k]
        bestAnswers.append(first)
    for j in tqdm(range(n)):
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
            #replace set of motifs if score improved when we repeat for range of n
            if(score(motifs) < score(bestMotifs)):
                bestMotifs = motifs
        #replace set of motifs if score improved at end of each of the 1000 ierations within the for loop 
        if(score(bestMotifs) < score(bestAnswers)):
            bestAnswers = bestMotifs;
    return(bestAnswers)

#function takes in a set of sequences and produces a matrix of rows based on the length of the sequences and 4 columns, each corresponding to the 4 nucleotides: A, C, G, and T 
#each column is a probability distribution of the likelihoods of each nucleotide appearing in that spot based on the motif matrix we are given.
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

#produces a random integer between 0 and n
def rand(n: int) -> int: 
    return(random.randint(0, n-1))

#takes in an input of a set of motifs and produces a score for that motif set, and we want to minimize the score
#score is assigned to each column based on the number of nucleotides in that column minus the most frequent nucleotide showing up in that column.
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

#Takes inputs of a (long) sequence, an integer k for k-mer size we want, and the profile of the motif set
#the code uses these inputs to go through all possible k-mers from the sequence, and calculate which is the most probable k-mer based on the profile matrix given
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

#Takes inputs of a string called pattern and the profile matrix based on a set of motifs and uses the probability distribution of each column to calculate the total probability of pattern occurring
#The function returns this probability and the function is used to find the profile most probable kmer.
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


    
