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
strainVCF = sys.argv[2]
outputName = sys.argv[3]


refDict = {}
strainDict = {}
chrListRef = []

#read in fasta
#build dictionary with chr as key and sequence as value

refGenome =  open(refGenomeName, 'r')
for seq_record in SeqIO.parse(refGenome, "fasta"):
    refDict[seq_record.id] = seq_record.seq
    idStr = str(seq_record.id)
    chrListRef.append(idStr)
    seqStr = str(seq_record.seq)
    strainDict[idStr] = MutableSeq(seqStr, IUPAC.IUPACAmbiguousDNA())
for k in strainDict:
    print(k)

outIndelName = outputName + ".indel"
outIndelFile = open(outIndelName, 'w')

vcf = open(strainVCF, 'r')
lines = vcf.readlines()
counter = 0
for line in lines:
    currentLine = line.strip('\n')
    if re.match('^#', currentLine): #Looks for lines that start with #
        header = currentLine
        #print "this is just the header"
    else:
        SNPline = currentLine.split() #Split based on tabs
        chr = SNPline[0]
        pos = int(SNPline[1])-1
        ref = SNPline[3]
        alt = SNPline[4]
        geno = SNPline[9]
        genoSplit = geno.split(":")
        if len(genoSplit) == 5: ##For some reason some are missing depth data so I'm just going to mask there
            genotype = genoSplit[0]
            alleleDepth = genoSplit[1]
            coverage = genoSplit[2]
            genoQual = genoSplit[3]
            alleleLikely = genoSplit[4]
            if coverage > 3:
                #print "Yay!"
                if len(ref) == 1 and len(alt) == 1: #Check that they arn't indels
                    if re.match('0', geno):
                        if (re.match('A', ref) and re.match('G', alt)) or (re.match('G', ref) and re.match('A', alt)):
                            alt = 'R'
                        if (re.match('C', ref) and re.match('T', alt)) or (re.match('T', ref) and re.match('C', alt)):
                            alt = 'Y'
                        if (re.match('G', ref) and re.match('T', alt)) or (re.match('G', ref) and re.match('T', alt)):
                            alt = 'K'
                        if (re.match('A', ref) and re.match('C', alt)) or (re.match('C', ref) and re.match('A', alt)):
                            alt = 'M'
                        if (re.match('C', ref) and re.match('G', alt)) or (re.match('G', ref) and re.match('C', alt)):
                            alt = 'S'
                        if (re.match('A', ref) and re.match('T', alt)) or (re.match('T', ref) and re.match('A', alt)):
                            alt = 'W'
                    else:
                        alt = alt
                    #print(chr+" "+str(pos)+" "+alt)
                    print(chr+"_"+str(len(strainDict[chr]))+"_"+str(pos))
                    sys.stdout.flush()
                    strainDict[chr][pos] = alt
                else:
                    #print "INDEL"
                    genotype =str()
                    if re.match('0', geno):
                        genotype = 'heterozygous'
                    else:
                        genotype = 'homozygous'
                    outIndelFile.write(chr+"\t"+str(pos)+"\t"+ref+"\t"+alt+"\t"+genotype+"\n")
                    replaceLen = len(ref)
                    if (replaceLen == 1):
                        replaceLen = len(alt)
                    if pos+replaceLen > len(strainDict[chr]):
                        strainDict[chr][pos:(len(strainDict[chr]))] = 'N'*(len(strainDict[chr])-pos)
                    else:
                        repN = 'N'*replaceLen
                        strainDict[chr][pos:(pos+replaceLen)] = repN
            else:
                strainDict[chr][pos] = 'N'
        else:
            strainDict[chr][pos] = 'N'
            counter = counter + 1


#print (str(counter)+" positions were missing depth data")
outputFASTAname = outputName + ".fasta"
outFile = open(outputFASTAname, "w")
#chrList = open(chrListname, 'r')
#chrLines = chrList.readlines()
for chr in chrListRef:
    #print(chr)
    outFile.write(">"+ chr + "\n")
    outFile.write(str(strainDict[chr]) + "\n")

#for line in chrLines:
#    chr = line.strip()
#    outFile.write(">"+ chr + "\n")
#    outFile.write(str(strainDict[chr]) + "\n")






