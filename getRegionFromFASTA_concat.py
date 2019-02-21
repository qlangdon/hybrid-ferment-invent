#!/usr/bin/env python
from Bio import SeqIO
import sys

fastaName = sys.argv[1]
genePosName = sys.argv[2]
outputPrefix = sys.argv[3]
outputName = outputPrefix+".fasta"
fastaFile = open(fastaName, 'r')
strainList = list()
strainDict = {}
genomeDict = {}
for seq_record in SeqIO.parse(fastaFile, "fasta"):
    strainList.append(seq_record.id)
    genomeDict[seq_record.id] = seq_record.seq
    strainDict[seq_record.id] = str()

with open(genePosName, 'r') as genePos:
    header = genePos.readline().rstrip('\n')
    for line in genePos.readlines():
        line = line.strip()
        # print(str(line))
        tokens = line.split()
        compStrain = tokens[0]
        startPos = int(tokens[1])
        endPos = int(tokens[2])
        for strain in strainList:
            region = str(genomeDict[strain][startPos + 1:endPos + 1])
            strainDict[strain] = strainDict[strain] + region

outFile = open(outputName, 'w')
for strain in strainList:
    outFile.write(">"+strain+"\n")
    outFile.write(strainDict[strain]+"\n")
outFile.close()
