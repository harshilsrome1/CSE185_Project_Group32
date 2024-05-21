import random 
from tqdm import tqdm

def gibbs_sampler(dna: list[str], k: int, t: int, n: int) -> list[str]:
    bestAnswers = []
    for eachdna in dna:
        number = rand(len(eachdna)-k)
        first = eachdna[number:number+k]
        bestAnswers.append(first)
    for j in tqdm(range(1000)):
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


    
