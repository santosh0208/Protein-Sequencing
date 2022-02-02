"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file= open(filename,"r")
    a=file.read()
    str=""
    for line in a.splitlines():
        str= str+line
    return str


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    condonlist = []
    var = ["UGA","UAG","UAA"]
    for word in range(startIndex,len(dna),3):
        dna = dna.replace("T","U")
        condon = dna[word:word+3]
        if condon not in var:
            condonlist.append(condon)
        else:
            condonlist.append(condon)
            break
    return condonlist 


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    File = open(filename,"r")
    protns = json.load(File)
    codondict= {}
    for key in protns:
        for values in protns[key]:
            values = values.replace("T","U")
            codondict[values] = key
    return codondict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    proteinlist = []
    for rna in codons:
        for rnaproteins in codonD:
            if rna == rnaproteins:
                proteinlist.append(codonD[rnaproteins])
                if proteinlist[0] == "Met":
                    proteinlist[0] = "Start"
    return proteinlist 


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    rnalist = readFile(dnaFilename)
    proteinlist = makeCodonDictionary(codonFilename)
    totalproteinlist = []
    j = 0
    unusedltrs = 0
    while j<len(rnalist) :
        word = rnalist[j:j+3]
        if word == "ATG":
            rna = dnaToRna(rnalist, j)
            totalproteinlist.append(generateProtein(rna,proteinlist))
            j = j+3*len(rna)
        else:
            unusedltrs += 1
            j += 1      
    return totalproteinlist


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    commonprotein = []
    for protein1 in proteinList1:
        for protein2 in proteinList2:
            if protein1 == protein2 and protein1 not in commonprotein:
                commonprotein.append(protein1)
    return commonprotein

'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    proteinlist = []
    for list in proteinList:
        for proteins in list:
            proteinlist.append(proteins)
    return proteinlist


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aminoaciddict = {}
    for word in aaList:
        if word not in aminoaciddict:
            aminoaciddict[word] = 0
        if word in aminoaciddict:
            aminoaciddict[word] += 1
    return aminoaciddict 


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    aminofreqlist =[]
    totalcount1 = len(combineProteins(proteinList1))
    totalcount2 = len(combineProteins(proteinList2))
    word = combineProteins(proteinList1) 
    word1 = combineProteins(proteinList2)
    count = aminoAcidDictionary(word)
    count1 = aminoAcidDictionary(word1)
    totalwrd = word + word1
    totalwrd = list(set(totalwrd))
    for i in totalwrd:
        freq1 = 0 
        freq2 = 0
        if i != "Start" and i != "Stop":
            if i not in aminofreqlist:
                if i in word and i in word1:
                    freq1 = count[i]/totalcount1
                    freq2 = count1[i]/totalcount2
                elif i in word and i not in word1:
                    freq1 = count[i]/totalcount1
                    freq2 = 0
                elif i not in word and i in word1:
                    freq1 = 0 
                    freq2 = count1[i]/totalcount2
                freqdiff = freq1-freq2
                if freqdiff > cutoff or freqdiff < -cutoff:
                    aminofreqlist.append([i,freq1,freq2])
    return aminofreqlist


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    list = []
    list1 = []
    str = ''
    for j in commonalities:
        if j not in list:
            list.append(j[1:-1])
    for k in list:
        if len(k) > 0:
            k = '-'.join(k)
            list1.append(k)
    list1.sort()
    for l in list1:
        str += ' '+ l +"\n"
    print(str)
    for b in differences:
        wrd = b[0]
        seq1 = round(b[1]*100,2)
        seq2 = round(b[2]*100,2)
        print(f"{wrd}:{seq1}% in Seq1, {seq2}% in seq2") 

    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    labels = []
    word = combineProteins(proteinList1) + combineProteins(proteinList2)
    for j in word:
        if j not in labels:
            labels.append(j)
    totallabels = sorted(labels)
    return totallabels


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    freqlist = []
    word = combineProteins(proteinList)
    totalcount = len(combineProteins(proteinList))
    count = aminoAcidDictionary(word)
    for j in labels:
        if j in count:
            freq = count[j]/totalcount
            freqlist.append(freq)
        else:
            freqlist.append(0)
    return freqlist


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    w = 0.35
    plt.bar(xLabels, freqList1, width=-w, align='edge', label=label1, edgecolor = edgeList)
    plt.bar(xLabels, freqList2, width= w, align='edge', label=label2, edgecolor = edgeList)

    plt.xticks(rotation="vertical")
    plt.legend()
    plt.title("Compare Human and Elephant Genes")

    plt.show()
    return 


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##

    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()


    ## Uncomment these for Week 3 ##

    # print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    # test.week3Tests()
    # print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    # runFullProgram()
