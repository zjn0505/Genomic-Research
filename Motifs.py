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
            profile[symbol].append(count[symbol][j] / t)
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
