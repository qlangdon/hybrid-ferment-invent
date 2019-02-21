#!/usr/bin/env python
import Bio
from Bio import SeqIO
import sys
import re
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

#Python script to extract data from fast and vcf
#print ('Arguments: refGenome, chromosomeList, strainVCF, output name  \n')

refGenomeName = sys.argv[1]
#chrListname = sys.argv[2]
phaseInfoFile = sys.argv[2]
outputName = sys.argv[3]


refDict = {}
hap1dict = {}
hap2dict = {}
chrListRef = []

#read in fasta
#build dictionary with chr as key and sequence as value

refGenome = open(refGenomeName, 'r')
for seq_record in SeqIO.parse(refGenome, "fasta"):
    refDict[seq_record.id] = seq_record.seq
    idStr = str(seq_record.id)
    chrListRef.append(idStr)
    seqStr = str(seq_record.seq)
    hap1dict[idStr] = MutableSeq(seqStr, IUPAC.IUPACAmbiguousDNA())
    hap2dict[idStr] = MutableSeq(seqStr, IUPAC.IUPACAmbiguousDNA())
#for k in strainDict:
#    print(k)
refGenome.close()

#counter = 0
phaseInfo = open(phaseInfoFile, 'r')
lines = phaseInfo.readlines()[1:]
for line in lines:
    currentLine = line.strip('\n')
    alleleInfo = currentLine.split('\t')
    chr = alleleInfo[0]
    chrPos = int(alleleInfo[1])-1
    depth = int(alleleInfo[3])
    het = alleleInfo[4]
    ref = alleleInfo[5]
    refType = alleleInfo[6]
    alt = alleleInfo[8]
    altType = alleleInfo[9]
    firstHap = alleleInfo[12]
    secondHap = alleleInfo[13]
    if depth > 3:
        if refType == "SNP" and altType == "SNP":
            if het == 'TRUE':
                if int(firstHap) == 1:
                    hap1nt = ref
                    hap2nt = alt
                elif int(firstHap) == 2:
                    hap1nt = alt
                    hap2nt = ref
                hap1dict[chr][chrPos] = hap1nt
                hap2dict[chr][chrPos] = hap2nt
            else:
                hap1dict[chr][chrPos] = alt
                hap2dict[chr][chrPos] = alt
        elif refType == "del" or altType == "in":
            replaceLen = len(ref)
            if replaceLen == 1:
                replaceLen = len(alt)
            if chrPos + replaceLen > len(hap1dict[chr]):
                hap1dict[chr][chrPos:(len(hap1dict[chr]))] = 'N' * (len(hap1dict[chr]) - chrPos)
                hap2dict[chr][chrPos:(len(hap2dict[chr]))] = 'N' * (len(hap2dict[chr]) - chrPos)
            else:
                repN = 'N' * replaceLen
                hap1dict[chr][chrPos:(chrPos + replaceLen)] = repN
                hap2dict[chr][chrPos:(chrPos + replaceLen)] = repN
    else:
        hap1dict[chr][chrPos] = 'N'
        hap2dict[chr][chrPos] = 'N'

phaseInfo.close()

#print (str(counter)+" positions were missing depth data")
hap1FASTAname = outputName + "_phased_1.fasta"
hap2FASTAname = outputName + "_phased_2.fasta"
hap1fasta = open(hap1FASTAname, "w")
hap2fasta = open(hap2FASTAname, "w")
for chr in chrListRef:
    #print(chr)
    hap1fasta.write(">"+ chr + "\n")
    hap2fasta.write(">" + chr + "\n")
    hap1fasta.write(str(hap1dict[chr]) + "\n")
    hap2fasta.write(str(hap2dict[chr]) + "\n")

hap1fasta.close()
hap2fasta.close()

#for line in chrLines:
#    chr = line.strip()
#    outFile.write(">"+ chr + "\n")
#    outFile.write(str(strainDict[chr]) + "\n")






