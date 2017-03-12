# Copy your updated FrequentWords function (along with all required subroutines) below this line
# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates 

# Input:  A list Items
# Output: A list containing all objects from Items without duplicates
def remove_duplicates(Items):
    ItemsNoDuplicates = [] # output variable
    # your code here
    for item in Items:
        if item not in ItemsNoDuplicates:
            ItemsNoDuplicates.append(item)
            
    return ItemsNoDuplicates

# Input:  A string Text and an integer k
# Output: CountDict(Text, k)
# HINT:   This code should be identical to when you last implemented CountDict
def CountDict(Text, k):
    Count = {} # output variable
    # your code here
    for i in range(len(Text) - k + 1):
        Pattern = Text[i : i + k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

# Input:  Strings Pattern and Text
# Output: The number of times Pattern appears in Text
# HINT:   This code should be identical to when you last implemented PatternCount
def PatternCount(Pattern, Text):
    count = 0 # output variable
    # your code here
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i : i + len(Pattern)] == Pattern:
            count += 1
    return count
	
# Input:  Two strings, Pattern and Genome
# Output: A list containing all starting positions where Pattern appears as a substring of Genome
def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    # your code here
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i : i + len(Pattern)] == Pattern:
            positions.append(i)
    return positions



# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    # your code here
    count = 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            count += 1
    return count


# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions = [] # output variable
    # your code here
    skew = Skew(Genome)
    minValue = min(skew.values())
    for i in range(0, len(skew)):
        if skew[i] == minValue:
            positions.append(i)
    return positions

# Input:  A String Genome
# Output: Skew(Genome)
# HINT:   This code should be taken from the last Code Challenge.
def Skew(Genome):
    skew = {} #initializing the dictionary
    # your code here
    skew[0] = 0
    for i in range(0, len(Genome)):
        if Genome[i] == "C":
            skew[i + 1] = skew[i] - 1
        elif Genome[i] == "G":
            skew[i + 1] = skew[i] + 1
        else:
            skew[i + 1] = skew[i]
    return skew

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # initializing list of positions
    # your code here
    for i in range(0, len(Text)-len(Pattern) + 1):
        chop = Text[i: i + len(Pattern)]
        if HammingDistance(Pattern, chop) <= d:
            positions.append(i)
    return positions

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    # your code here
    for i in range(0, len(Text)-len(Pattern) + 1):
        chop = Text[i: i + len(Pattern)]
        if HammingDistance(Pattern, chop) <= d:
            count += 1
    return count

# Now set Text equal to the Vibrio cholerae oriC and k equal to 10

oriC = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"

k = 10
# Finally, print the result of calling FrequentWords on Text and k.

print(FrequentWords(oriC, k))