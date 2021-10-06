import collections
from Bio import SeqIO
import numpy as np

aminoacidDict = {
    'ttt': 'Phe', 'tct': 'Ser', 'tat': 'Tyr', 'tgt': 'Cys',
    'ttc': 'Phe', 'tcc': 'Ser', 'tac': 'Tyr', 'tgc': 'Cys',
    'tta': 'Leu', 'tca': 'Ser', 'taa': '*', 'tga': '*',
    'ttg': 'Leu', 'tcg': 'Ser', 'tag': '*', 'tgg': 'Trp',
    'ctt': 'Leu', 'cct': 'Pro', 'cat': 'His', 'cgt': 'Arg',
    'ctc': 'Leu', 'ccc': 'Pro', 'cac': 'His', 'cgc': 'Arg',
    'cta': 'Leu', 'cca': 'Pro', 'caa': 'Gln', 'cga': 'Arg',
    'ctg': 'Leu', 'ccg': 'Pro', 'cag': 'Gln', 'cgg': 'Arg',
    'att': 'Ile', 'act': 'Thr', 'aat': 'Asn', 'agt': 'Ser',
    'atc': 'Ile', 'acc': 'Thr', 'aac': 'Asn', 'agc': 'Ser',
    'ata': 'Ile', 'aca': 'Thr', 'aaa': 'Lys', 'aga': 'Arg',
    'atg': 'Met', 'acg': 'Thr', 'aag': 'Lys', 'agg': 'Arg',
    'gtt': 'Val', 'gct': 'Ala', 'gat': 'Asp', 'ggt': 'Gly',
    'gtc': 'Val', 'gcc': 'Ala', 'gac': 'Asp', 'ggc': 'Gly',
    'gta': 'Val', 'gca': 'Ala', 'gaa': 'Glu', 'gga': 'Gly',
    'gtg': 'Val', 'gcg': 'Ala', 'gag': 'Glu', 'ggg': 'Gly'
}

listOfCodons = {
    'ttt', 'ttc', 'tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg', 'att', 'atc', 'ata', 'atg', 'gtt',
    'gtc', 'gta', 'gtg', 'tct', 'tcc', 'tca', 'tcg', 'cct', 'ccc', 'cca', 'ccg', 'act', 'acc',
    'aca', 'acg', 'gct', 'gcc', 'gca', 'gcg', 'tat', 'tac', 'taa', 'tag', 'cat', 'cac', 'caa',
    'cag', 'aat', 'aac', 'aaa', 'aag', 'gat', 'gac', 'gaa', 'gag', 'tgt', 'tgc', 'tga',
    'tgg', 'cgt', 'cgc', 'cga', 'cgg', 'agt', 'agc', 'aga', 'agg', 'ggt', 'ggc', 'gga', 'ggg'
}


def makingListOfDicodons():
    listOfDicodons = []
    for a in listOfCodons:
        for b in listOfCodons:
            c = (a + b)
            listOfDicodons.append(c)
    return listOfDicodons


listOfDicodons = makingListOfDicodons()
listOfDicodons = dict.fromkeys(listOfDicodons, 0)

listOfUniqueCodons = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Tyr', '*', 'His', 'Gln', 'Asn',
                      'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly']


# nuskaito failą
def readFile(file_name):
    for sequence_record in SeqIO.parse(open("data/" + file_name), "fasta"):
        return sequence_record


# randa kodonus su ATG pradžia ir su TAA, TAG ar TGA pabaiga
def findCodonSequence(sequence):
    i = 0
    codons_list = []
    while i < len(sequence):
        if sequence[i] == 'ATG':
            start_position = i
            j = i
            while j < len(sequence):
                if sequence[j] == 'TAA' or sequence[j] == 'TAG' or sequence[j] == 'TGA':
                    end_position = j
                    codon = ''.join(str(triplet) for triplet in sequence[start_position: end_position + 1])
                    codons_list.append(codon)
                    i = j
                    break
                j = j + 1
        i = i + 1
    codons_list = filterCodons(codons_list)  # atfiltruoja sekos
    return codons_list


# suskaido seką į tripletus ir frames bei po to padaro vieną didelę seką iš kodonų
def getTriplets(sequence_record):
    listFrame1 = [sequence_record.seq[i:i + 3] for i in range(0, len(sequence_record.seq), 3)]
    listFrame2 = [sequence_record.seq[i:i + 3] for i in range(1, len(sequence_record.seq), 3)]
    listFrame3 = [sequence_record.seq[i:i + 3] for i in range(2, len(sequence_record.seq), 3)]
    listReverseFrame1 = [sequence_record.seq.reverse_complement()[i:i + 3] for i in
                         range(0, len(sequence_record.seq), 3)]
    listReverseFrame2 = [sequence_record.seq.reverse_complement()[i:i + 3] for i in
                         range(1, len(sequence_record.seq), 3)]
    listReverseFrame3 = [sequence_record.seq.reverse_complement()[i:i + 3] for i in
                         range(2, len(sequence_record.seq), 3)]

    return findCodonSequence(listFrame1) + findCodonSequence(listFrame2) + findCodonSequence(
        listFrame3) + findCodonSequence(listReverseFrame1) + findCodonSequence(listReverseFrame2) + findCodonSequence(
        listReverseFrame3)


# suranda ilgiausią kodoną (ilgiausias atstumas nuo stop kodono iki start kodono)
def getLongestCodon(codons_list):
    longestCodon = max(codons_list, key=len)
    return longestCodon


# atfiltruoja visus fragmentus ir palieka tik tuos, kurie ilgesni negu 100
def filterCodons(codons_list):
    filtered = [codon for codon in codons_list if len(codon) >= 100]
    return filtered


# randa kodono dažnį
def findCodonFrequency(sequence):
    codon_data = []
    sequenceString = ''.join(sequence)
    sequenceList = sequenceString.lower()
    sequenceTriplets = [sequenceList[i:i + 3] for i in range(0, len(sequenceList), 3)]
    codon_count = collections.Counter(aminoacidDict[c] for c in sequenceTriplets)
    codon_number = len(sequenceTriplets)
    for codon in codon_count:
        codon_data = np.append(codon_data, [round(codon_count[codon] / codon_number, 4)], axis=0)
    return codon_data


# randa dikodono dažnį
def findDicodonFrequency(sequence):
    dicodon_data = []
    sequenceString = ''.join(sequence)
    sequenceList = sequenceString.lower()
    dicodons = [str(sequenceList[i:i + 6]) for i in range(0, len(sequenceList) - 3, 3)]
    dicodon_count = collections.Counter(dicodons)
    updateDict = listOfDicodons
    updateDict.update(dicodon_count)
    frequency = list(updateDict.values())
    dicodon_number = len(dicodons)
    for dicodon in frequency:
        dicodon_data = np.append(dicodon_data, round(frequency[dicodon]/dicodon_number, 4))
    return dicodon_data


n = 64
m = 4096

bacterial1 = readFile("bacterial1.fasta")
bacterial2 = readFile("bacterial2.fasta")
bacterial3 = readFile("bacterial3.fasta")
bacterial4 = readFile("bacterial4.fasta")
mamalian1 = readFile("mamalian1.fasta")
mamalian2 = readFile("mamalian2.fasta")
mamalian3 = readFile("mamalian3.fasta")
mamalian4 = readFile("mamalian4.fasta")

bacterial1_codon = findCodonFrequency(getTriplets(bacterial1))
bacterial2_codon = findCodonFrequency(getTriplets(bacterial2))
bacterial3_codon = findCodonFrequency(getTriplets(bacterial3))
bacterial4_codon = findCodonFrequency(getTriplets(bacterial4))
mamalian1_codon = findCodonFrequency(getTriplets(mamalian1))
mamalian2_codon = findCodonFrequency(getTriplets(mamalian2))
mamalian3_codon = findCodonFrequency(getTriplets(mamalian3))
mamalian4_codon = findCodonFrequency(getTriplets(mamalian4))

bacterial1_dicodon = findDicodonFrequency(getTriplets(bacterial1))
bacterial2_dicodon = findDicodonFrequency(getTriplets(bacterial2))
bacterial3_dicodon = findDicodonFrequency(getTriplets(bacterial3))
bacterial4_dicodon = findDicodonFrequency(getTriplets(bacterial4))
mamalian1_dicodon = findDicodonFrequency(getTriplets(mamalian1))
mamalian2_dicodon = findDicodonFrequency(getTriplets(mamalian2))
mamalian3_dicodon = findDicodonFrequency(getTriplets(mamalian3))
mamalian4_dicodon = findDicodonFrequency(getTriplets(mamalian4))

bacterial1_codon_distance = []
bacterial2_codon_distance = []
bacterial3_codon_distance = []
bacterial4_codon_distance = []
mamalian1_codon_distance = []
mamalian2_codon_distance = []
mamalian3_codon_distance = []
mamalian4_codon_distance = []

bacterial1_dicodon_distance = []
bacterial2_dicodon_distance = []
bacterial3_dicodon_distance = []
bacterial4_dicodon_distance = []
mamalian1_dicodon_distance = []
mamalian2_dicodon_distance = []
mamalian3_dicodon_distance = []
mamalian4_dicodon_distance = []

#kodonu matricų sudarimas
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - bacterial1_codon), 2)) / n))
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - bacterial2_codon), 2)) / n))
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - bacterial3_codon), 2)) / n))
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - bacterial4_codon), 2)) / n))
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - mamalian1_codon), 2)) / n))
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - mamalian2_codon), 2)) / n))
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - mamalian3_codon), 2)) / n))
bacterial1_codon_distance.append(np.sqrt(np.sum(np.power((bacterial1_codon - mamalian4_codon), 2)) / n))

bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - bacterial1_codon), 2)) / n))
bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - bacterial2_codon), 2)) / n))
bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - bacterial3_codon), 2)) / n))
bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - bacterial4_codon), 2)) / n))
bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - mamalian1_codon), 2)) / n))
bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - mamalian2_codon), 2)) / n))
bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - mamalian3_codon), 2)) / n))
bacterial2_codon_distance.append(np.sqrt(np.sum(np.power((bacterial2_codon - mamalian4_codon), 2)) / n))

bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - bacterial1_codon), 2)) / n))
bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - bacterial2_codon), 2)) / n))
bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - bacterial3_codon), 2)) / n))
bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - bacterial4_codon), 2)) / n))
bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - mamalian1_codon), 2)) / n))
bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - mamalian2_codon), 2)) / n))
bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - mamalian3_codon), 2)) / n))
bacterial3_codon_distance.append(np.sqrt(np.sum(np.power((bacterial3_codon - mamalian4_codon), 2)) / n))

bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - bacterial1_codon), 2)) / n))
bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - bacterial2_codon), 2)) / n))
bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - bacterial3_codon), 2)) / n))
bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - bacterial4_codon), 2)) / n))
bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - mamalian1_codon), 2)) / n))
bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - mamalian2_codon), 2)) / n))
bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - mamalian3_codon), 2)) / n))
bacterial4_codon_distance.append(np.sqrt(np.sum(np.power((bacterial4_codon - mamalian4_codon), 2)) / n))

mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - bacterial1_codon), 2)) / n))
mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - bacterial2_codon), 2)) / n))
mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - bacterial3_codon), 2)) / n))
mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - bacterial4_codon), 2)) / n))
mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - mamalian1_codon), 2)) / n))
mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - mamalian2_codon), 2)) / n))
mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - mamalian3_codon), 2)) / n))
mamalian1_codon_distance.append(np.sqrt(np.sum(np.power((mamalian1_codon - mamalian4_codon), 2)) / n))

mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - bacterial1_codon), 2)) / n))
mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - bacterial2_codon), 2)) / n))
mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - bacterial3_codon), 2)) / n))
mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - bacterial4_codon), 2)) / n))
mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - mamalian1_codon), 2)) / n))
mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - mamalian2_codon), 2)) / n))
mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - mamalian3_codon), 2)) / n))
mamalian2_codon_distance.append(np.sqrt(np.sum(np.power((mamalian2_codon - mamalian4_codon), 2)) / n))

mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - bacterial1_codon), 2)) / n))
mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - bacterial2_codon), 2)) / n))
mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - bacterial3_codon), 2)) / n))
mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - bacterial4_codon), 2)) / n))
mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - mamalian1_codon), 2)) / n))
mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - mamalian2_codon), 2)) / n))
mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - mamalian3_codon), 2)) / n))
mamalian3_codon_distance.append(np.sqrt(np.sum(np.power((mamalian3_codon - mamalian4_codon), 2)) / n))

mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - bacterial1_codon), 2)) / n))
mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - bacterial2_codon), 2)) / n))
mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - bacterial3_codon), 2)) / n))
mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - bacterial4_codon), 2)) / n))
mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - mamalian1_codon), 2)) / n))
mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - mamalian2_codon), 2)) / n))
mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - mamalian3_codon), 2)) / n))
mamalian4_codon_distance.append(np.sqrt(np.sum(np.power((mamalian4_codon - mamalian4_codon), 2)) / n))

#dikodonų matricų sudarymas
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - bacterial1_dicodon), 2)) / m))
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - bacterial2_dicodon), 2)) / m))
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - bacterial3_dicodon), 2)) / m))
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - bacterial4_dicodon), 2)) / m))
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - mamalian1_dicodon), 2)) / m))
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - mamalian2_dicodon), 2)) / m))
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - mamalian3_dicodon), 2)) / m))
bacterial1_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial1_dicodon - mamalian4_dicodon), 2)) / m))

bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - bacterial1_dicodon), 2)) / m))
bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - bacterial2_dicodon), 2)) / m))
bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - bacterial3_dicodon), 2)) / m))
bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - bacterial4_dicodon), 2)) / m))
bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - mamalian1_dicodon), 2)) / m))
bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - mamalian2_dicodon), 2)) / m))
bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - mamalian3_dicodon), 2)) / m))
bacterial2_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial2_dicodon - mamalian4_dicodon), 2)) / m))

bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - bacterial1_dicodon), 2)) / m))
bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - bacterial2_dicodon), 2)) / m))
bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - bacterial3_dicodon), 2)) / m))
bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - bacterial4_dicodon), 2)) / m))
bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - mamalian1_dicodon), 2)) / m))
bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - mamalian2_dicodon), 2)) / m))
bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - mamalian3_dicodon), 2)) / m))
bacterial3_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial3_dicodon - mamalian4_dicodon), 2)) / m))

bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - bacterial1_dicodon), 2)) / m))
bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - bacterial2_dicodon), 2)) / m))
bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - bacterial3_dicodon), 2)) / m))
bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - bacterial4_dicodon), 2)) / m))
bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - mamalian1_dicodon), 2)) / m))
bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - mamalian2_dicodon), 2)) / m))
bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - mamalian3_dicodon), 2)) / m))
bacterial4_dicodon_distance.append(np.sqrt(np.sum(np.power((bacterial4_dicodon - mamalian4_dicodon), 2)) / m))

mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - bacterial1_dicodon), 2)) / m))
mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - bacterial2_dicodon), 2)) / m))
mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - bacterial3_dicodon), 2)) / m))
mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - bacterial4_dicodon), 2)) / m))
mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - mamalian1_dicodon), 2)) / m))
mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - mamalian2_dicodon), 2)) / m))
mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - mamalian3_dicodon), 2)) / m))
mamalian1_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian1_dicodon - mamalian4_dicodon), 2)) / m))

mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - bacterial1_dicodon), 2)) / m))
mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - bacterial2_dicodon), 2)) / m))
mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - bacterial3_dicodon), 2)) / m))
mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - bacterial4_dicodon), 2)) / m))
mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - mamalian1_dicodon), 2)) / m))
mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - mamalian2_dicodon), 2)) / m))
mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - mamalian3_dicodon), 2)) / m))
mamalian2_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian2_dicodon - mamalian4_dicodon), 2)) / m))

mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - bacterial1_dicodon), 2)) / m))
mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - bacterial2_dicodon), 2)) / m))
mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - bacterial3_dicodon), 2)) / m))
mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - bacterial4_dicodon), 2)) / m))
mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - mamalian1_dicodon), 2)) / m))
mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - mamalian2_dicodon), 2)) / m))
mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - mamalian3_dicodon), 2)) / m))
mamalian3_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian3_dicodon - mamalian4_dicodon), 2)) / m))

mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - bacterial1_dicodon), 2)) / m))
mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - bacterial2_dicodon), 2)) / m))
mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - bacterial3_dicodon), 2)) / m))
mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - bacterial4_dicodon), 2)) / m))
mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - mamalian1_dicodon), 2)) / m))
mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - mamalian2_dicodon), 2)) / m))
mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - mamalian3_dicodon), 2)) / m))
mamalian4_dicodon_distance.append(np.sqrt(np.sum(np.power((mamalian4_dicodon - mamalian4_dicodon), 2)) / m))

file = open("codon_distance_matrix.txt", "w")
file.write('8\n')
file.write('Lactococcus '+str(bacterial1_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Escherichia '+str(bacterial2_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Streptococcus '+str(bacterial3_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Cellulophaga '+str(bacterial4_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Coronavirus '+str(mamalian1_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Adenovirus '+str(mamalian2_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Variolavirus '+str(mamalian3_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Herpesvirus '+str(mamalian4_codon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.close()

file = open("dicodon_distance_matrix.txt", "w")
file.write('8\n')
file.write('Lactococcus '+str(bacterial1_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Escherichia '+str(bacterial2_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Streptococcus '+str(bacterial3_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Cellulophaga '+str(bacterial4_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Coronavirus '+str(mamalian1_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Adenovirus '+str(mamalian2_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Variolavirus '+str(mamalian3_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.write('Herpesvirus '+str(mamalian4_dicodon_distance).translate(dict.fromkeys(map(ord, u"[],")))+'\n')
file.close()
