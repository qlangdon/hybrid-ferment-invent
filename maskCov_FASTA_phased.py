__author__ = 'Quinn'

#!/usr/bin/env python
from Bio import SeqIO
import sys
import re
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

#Python script to extract data from fast and vcf
#print 'Arguments: coverage quant by strain  \n'

strainCovQuantName = sys.argv[1]

maskingInfo = ""
strainCovQuant = open(strainCovQuantName, 'r')
strains = strainCovQuant.readlines()
strain = strains[1]
lowerCov = int(10)
upperCov = int(500)
strainLine = strain.strip('\n')
covInfo = strainLine.split()
strainName = re.sub(r'"','',covInfo[0])
hap1name = strainName + "_phased_1"
hap2name = strainName + "_phased_2"
lowerCheck = int(covInfo[1])
upperCheck = int(covInfo[2])
if lowerCheck < lowerCov:
    #print('Lower < 10')
    if lowerCheck < 3:
        lowerCov = int(3)
    else:
        lowerCov = lowerCheck
if upperCheck < lowerCov:
    upperCov = int(10)
else:
    upperCov = upperCheck
maskingInfo = maskingInfo + hap1name + ":\tlowLimit=" + str(lowerCov) + "\tupperLimit=" + str(upperCov) + "\n"

fastaToChangeName = hap1name+".fasta"
fastaToChange = open(fastaToChangeName, 'r')

strainDict = {}
genomeDict = {}
chrLength = {}
chrListRef = []
totalLen = float()
for seq_record in SeqIO.parse(fastaToChange, "fasta"):
    strainDict[seq_record.id] = seq_record.seq
    idStr = str(seq_record.id)
    seqStr = str(seq_record.seq)
    chrListRef.append(idStr)
    chrLength[idStr] = len(seqStr)
    totalLen = totalLen + float(len(seqStr))
    #print(len(seqStr))
    genomeDict[idStr] = MutableSeq(seqStr, IUPAC.IUPACAmbiguousDNA())

Ncount = 0

for key in genomeDict:
    Ncount = Ncount + genomeDict[key].count("N")

lenToMaskUpper=0
lenToMaskLower=0
existingN = 0
mito = 0
strainBed = strainName +".bedgraph"
bed = open(strainBed, 'r')
lines = bed.readlines()
for line in lines:
    currentLine = line.strip('\n')
    info = currentLine.split()
    chrom = info[0]
    chromstart = int(info[1])-1
    value = int(info[2])
    if value>upperCov:
        if re.match('Scer_chrmt', chrom):
            mito = mito +1
            #genomeLen = 23704987.0
        else:
            toReplace = genomeDict[chrom][chromstart]
            if re.match('N', toReplace):
                existingN = existingN + 1
            else:
                genomeDict[chrom][chromstart] = 'N'
        lenToMaskUpper = lenToMaskUpper + 1
    elif value<=lowerCov:
        if re.match('Scer_chrmt', chrom):
            mito = mito +1
        else:
            toReplace = genomeDict[chrom][chromstart]
            if re.match('N', toReplace):
                existingN = existingN + 1
            else:
                genomeDict[chrom][chromstart] = 'N'
        lenToMaskLower = lenToMaskLower + 1
totalCovMasked = lenToMaskLower + lenToMaskUpper
totalMasked = totalCovMasked + Ncount - existingN
rawNumMasked =  "\t\tExisting N=" + str(existingN) + "\tLower to Mask=" + str(lenToMaskLower) + "\tUpper to Mask=" + str(lenToMaskUpper)
percentMasked =  "\t\tExisting N=" + "%.2f" % ((existingN/totalLen)*100) + "%\tLower to Mask=" + "%.2f" % ((lenToMaskLower/totalLen)*100) + "%\tUpper to Mask=" + "%.2f" % ((lenToMaskUpper/totalLen)*100) + "%"
strTotalMasked = "\t\tRaw Total Masked=" + str(totalMasked) + "\tPercent Total Masked=" + "%.2f" % ((totalMasked/totalLen)*100) + "%"
#print(rawNumMasked)
#print(percentMasked)
#print(strTotalMasked)
maskingInfo = maskingInfo + rawNumMasked + "\n" + percentMasked + "\n" + strTotalMasked + "\n"

outputFileName = hap1name+"_covMasked.fasta"
outFile = open(outputFileName, 'w')
#if re.match('.+lager', strainName):
#    chrListName = "chromosome_combo.txt"
#else:
#    chrListName = "chromosome.txt"
#chrList = open(chrListName, 'r')
#chrLines = chrList.readlines()
for chr in chrListRef:
    outFile.write(">"+ chr + "\n")
    outFile.write(str(genomeDict[chr]) + "\n")
outFile.close()

maskingInfo = maskingInfo + hap2name + ":\tlowLimit=" + str(lowerCov) + "\tupperLimit=" + str(upperCov) + "\n"

fastaToChangeName = hap2name+".fasta"
fastaToChange = open(fastaToChangeName, 'r')

strainDict = {}
genomeDict = {}
chrLength = {}
chrListRef = []
totalLen = float()
for seq_record in SeqIO.parse(fastaToChange, "fasta"):
    strainDict[seq_record.id] = seq_record.seq
    idStr = str(seq_record.id)
    seqStr = str(seq_record.seq)
    chrListRef.append(idStr)
    chrLength[idStr] = len(seqStr)
    totalLen = totalLen + float(len(seqStr))
    #print(len(seqStr))
    genomeDict[idStr] = MutableSeq(seqStr, IUPAC.IUPACAmbiguousDNA())

Ncount = 0

for key in genomeDict:
    Ncount = Ncount + genomeDict[key].count("N")

lenToMaskUpper=0
lenToMaskLower=0
existingN = 0
mito = 0
strainBed = strainName +".bedgraph"
bed = open(strainBed, 'r')
lines = bed.readlines()
for line in lines:
    currentLine = line.strip('\n')
    info = currentLine.split()
    chrom = info[0]
    chromstart = int(info[1])-1
    value = int(info[2])
    if value>upperCov:
        if re.match('Scer_chrmt', chrom):
            mito = mito +1
            #genomeLen = 23704987.0
        else:
            toReplace = genomeDict[chrom][chromstart]
            if re.match('N', toReplace):
                existingN = existingN + 1
            else:
                genomeDict[chrom][chromstart] = 'N'
        lenToMaskUpper = lenToMaskUpper + 1
    elif value<=lowerCov:
        if re.match('Scer_chrmt', chrom):
            mito = mito +1
        else:
            toReplace = genomeDict[chrom][chromstart]
            if re.match('N', toReplace):
                existingN = existingN + 1
            else:
                genomeDict[chrom][chromstart] = 'N'
        lenToMaskLower = lenToMaskLower + 1
totalCovMasked = lenToMaskLower + lenToMaskUpper
totalMasked = totalCovMasked + Ncount - existingN
rawNumMasked =  "\t\tExisting N=" + str(existingN) + "\tLower to Mask=" + str(lenToMaskLower) + "\tUpper to Mask=" + str(lenToMaskUpper)
percentMasked =  "\t\tExisting N=" + "%.2f" % ((existingN/totalLen)*100) + "%\tLower to Mask=" + "%.2f" % ((lenToMaskLower/totalLen)*100) + "%\tUpper to Mask=" + "%.2f" % ((lenToMaskUpper/totalLen)*100) + "%"
strTotalMasked = "\t\tRaw Total Masked=" + str(totalMasked) + "\tPercent Total Masked=" + "%.2f" % ((totalMasked/totalLen)*100) + "%"
#print(rawNumMasked)
#print(percentMasked)
#print(strTotalMasked)
maskingInfo = maskingInfo + rawNumMasked + "\n" + percentMasked + "\n" + strTotalMasked + "\n"

outputFileName = hap2name+"_covMasked.fasta"
outFile = open(outputFileName, 'w')
#if re.match('.+lager', hap2name):
#    chrListName = "chromosome_combo.txt"
#else:
#    chrListName = "chromosome.txt"
#chrList = open(chrListName, 'r')
#chrLines = chrList.readlines()
for chr in chrListRef:
    outFile.write(">"+ chr + "\n")
    outFile.write(str(genomeDict[chr]) + "\n")
outFile.close()

infoFileName = strainName+"_phased_covMaskedInfo.txt"
infoFile = open(infoFileName, 'w')
infoFile.write(maskingInfo)
infoFile.close()



