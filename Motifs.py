# Count the occurence of 'ACGT' on each site
# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Add base 1 instead of 0 to every count, to prevent 0 in probililty.
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {}
    t = len(Motifs)
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Profile is the matrix for each site's possibility on ACGT
# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / float(t))
    return profile

# A better profile without 0, because there is chance that the sun may not raise tomorrow.
# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} # output variable
    count = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / float(t + 4))
    return profile

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
    
# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    count = Count(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    k = len(consensus)
    for i in range(k):
        for symbol in "ACGT":
            if symbol != consensus[i]:
                score += count[symbol][i]
    return score

# The Probility of input sequence given current profile
# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    probability = 1
    for i in range(len(Text)):
        probability *= Profile[Text[i]][i]
    return probability

# Find the most possible k-mer in Text given Profile    
# Input:  String Text, an integer k, and profile matrix Profile
# Output: ProfileMostProbablePattern(Text, k, Profile)
def ProfileMostProbablePattern(Text, k, Profile):
    best_prob = 0
    best_pattern = ""
    for i in range(len(Text) - k + 1):        
        s = Text[i:i+k]
        prob = Pr(s, Profile)
        if prob > best_prob:
            best_prob = prob
            best_pattern = s
    if best_prob == 0:
        best_pattern = Text[:k]
    return best_pattern

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Given Profile, find the Motif matrix from Dna
# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    results = []
    profileLen = len(Profile['A'])
    for dnaEntity in Dna:
        result = ProfileMostProbablePattern(dnaEntity, profileLen, Profile)
        results.append(result)
    return results

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
def RandomMotifs(Dna, k, t):
    results = []
    for i in range(0, t):
        rnd = random.randint(0, len(Dna[0])-k)
        result = Dna[i][rnd:rnd+k]
        results.append(result)
    return results

# Loop through a random motif to produce best motif
# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

# Repeat 1000 times to find the lowest score
def RepeatedRandomizedMotifSearch(Dna, k, t):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(1000):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs



# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    sum = 0
    for item in Probabilities:
        sum += Probabilities[item]
    for item in Probabilities:
        Probabilities[item] = Probabilities[item] / sum
    return Probabilities


# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    rnd = random.uniform(0,1)
    for item in list(Probabilities):
        rnd = rnd - Probabilities[item]
        if rnd < 0:
            return item

# Generate biased random substring with length k from Text given profile
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    Motifs = []
    for i in range(t):
        pos = random.randint(0, len(Dna[0])-k)
        Motifs.append(Dna[i][pos:pos+k])
    BestMotifs = Motifs
    for j in range(N):
        pick = random.randint(0, t-1)
        profile = ProfileWithPseudocounts(Motifs[0:i]+Motifs[i+1:t])
        Motifs[pick] = ProfileGeneratedString(Dna[pick], profile, k)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def RepeatedGibbsSampler(Dna, k, t, N):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(100):
        Motifs = GibbsSampler(Dna, k, t, N)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs
