__author__ = 'Quinn'

import sys
import re #to do regular expressions


vcfName = sys.argv[1] #vcf file name
outName = sys.argv[2]

vcfArray = list()
SNParray = list()
SNPinfo = list()
finalArray = list()
genotypeArray = list()
chrLen = dict()
chrList = list()
chrLenCumu = dict()
#totalHetCounter = 0
#totalHomAltCounter = 0
#homozygouteRefCounter = 0
#qualHetCounter = 0
#qualHomAltCounter = 0
#covHetCounter = 0
#covHomAltCounter = 0


outFile = open(outName+"_alleleInfo.txt", 'w')
outFile.write("chr\tchrPos\tgenomePos\tfilteredDepth\thetCalled\trefAllele\trefType\trefFreq\tallele_1\ttype_1\tfreq_1\tallele_2\ttype_2\tfreq_2\tallele_3\ttype_3\tfreq_3\tallele_4\ttype_4\tfreq_4\n")

vcf = open(vcfName, 'r')
lines = vcf.readlines()
for line in lines:
    currentLine = line.strip('\n')
    if re.match('##contig', currentLine): #look of the line with chr length info
        equalSplit = currentLine.split("=")
        commaSplit = equalSplit[2].split(",")
        chrName = commaSplit[0].strip("\s")
        chrLenStr = int(equalSplit[3].strip(">"))
        chrLen[chrName] = chrLenStr
        chrList.append(chrName)
        #SNParray.append(strainList)
    elif re.match('^#CHROM', currentLine):
        counter=0
        culumative = 0
        #print(chrList)
        for chr in chrList:
            chrLenCumu[chr] = culumative
            culumative = culumative+int(chrLen[chr])
    elif re.match('^\w', currentLine): #Look for the data lines, needs to be changed for other line starters
        SNPline = currentLine.split()  # Split based on tabs
        chr = SNPline[0]
        #print(chr)
        chrPos = SNPline[1]
        genomePos = chrLenCumu[chr] + int(SNPline[1])
        refAllele = SNPline[3]
        #print(refAllele)
        if len(refAllele) == 1:
            refType = "SNP"
        else:
            refType = "del"
        numAlt = 1
        altAllele = SNPline[4]
        altAlDict = {1:"NA",2:"NA",3:"NA",4:"NA"}
        altTyDict = {1:"NA",2:"NA",3:"NA",4:"NA"}
        if re.search(",",altAllele):
            altAllele = altAllele.split(",")
            numAlt = len(altAllele)
            allelNum = 1
            for i in range(0,len(altAllele)):
                altAlDict[allelNum] = altAllele[i]
                if len(altAllele[i])==1:
                    altTyDict[allelNum] = "SNP"
                else:
                    altTyDict[allelNum] = "in"
                allelNum = allelNum+1
        else:
            altAlDict[1] = altAllele
            if len(altAllele)==1:
                altTyDict[1] = "SNP"
            else:
                altTyDict[1] = "in"
        genotypeInfo = SNPline[9].split(":")
        if len(genotypeInfo) > 3:
            hetState = 'NA'
            stateCall = genotypeInfo[0]
            if stateCall=="0/1":
                hetState = 'TRUE'
            elif stateCall=="1/1":
                hetState='FALSE'
            DP = genotypeInfo[2]
            # print(DP)
            if int(DP) > 3:
                AD = genotypeInfo[1].split(",")
                refD = AD[0]
                #print(refD)
                refFreq = str(round((float(refD)/float(DP)),2))
                altD = AD[1:len(AD)]
                altDeDict = {1:"NA", 2:"NA", 3:"NA", 4:"NA"}
                if len(altD) == 1:
                    altDeDict[1] = str(round((float(altD[0]) / float(DP)), 2))
                else:
                    for i in range(0,len(altAllele)):
                        altDeDict[i+1] = str(round((float(altD[i])/float(DP)),2))
                outFile.write(chr + "\t" + chrPos + "\t" + str(genomePos) + "\t" + str(DP) + "\t" + hetState + "\t" + refAllele + "\t" + refType + "\t" + refFreq)
                for i in range(1,5,1):
                    outFile.write("\t" + altAlDict[i] + "\t" + altTyDict[i] + "\t" + altDeDict[i])
                outFile.write("\n")







#         counter = 9     #9th position is where genotype data starts
#         if len(SNPline[3]) == 1 and len(SNPline[4]) == 1: #Check that they aren't indels
#             SNPinfo = list()
#             SNPid = SNPline[0]+":"+SNPline[1]  #concatnate chr and position for SNP name
#             #print(SNPid)
#             chr = SNPline[0]
#             chrPos = SNPline[1]
#             genomePos = chrLenCumu[chr]+int(SNPline[1])
#             SNPinfo.append(SNPid)   #add it to the list
#             genotypeInfo = SNPline[9].split(":")
#             if re.match('\.', genotypeInfo[0]): #build the table with the SNP info
#                 SNPinfo.append(1) #if homozygous for the ref
#                 #outFile.write("1\t"+genotypeInfo[1]+"\t"+genotypeInfo[2])
#                 #print(counter - 9)
#                 #print(strainList[counter-9])
#                 homozygouteRefCounter = homozygouteRefCounter+1
#             elif re.match('0', genotypeInfo[0]):
#                 SNPinfo.append(2) #if heterozygous
#                 totalHetCounter = totalHetCounter+1
#                 if float(SNPline[5]) >= 30:
#                     qualHetCounter = qualHetCounter + 1
#                     if len(genotypeInfo)>3:
#                         if int(genotypeInfo[2]) >= 3:
#                             outFile.write(chr + "\t")
#                             outFile.write(chrPos + "\t")
#                             outFile.write(str(genomePos) + "\t")
#                             outFile.write("2\t")
#                             outFile.write("\n")
#                             if int(genotypeInfo[2]) >= 10:
#                                 covHetCounter = covHetCounter + 1
#                             #outFile.write(genotypeInfo[1]+"\t"+genotypeInfo[2])
#             elif re.match('^1', genotypeInfo[0]):
#                 SNPinfo.append(3) #if homozygous for the alternate allele
#                 totalHomAltCounter = totalHomAltCounter+1
#                 if float(SNPline[5]) >= 30:
#                     qualHomAltCounter = qualHomAltCounter + 1
#                     #print(genotypeInfo[2])
#                     if len(genotypeInfo)>3:
#                         if int(genotypeInfo[2]) >= 3:
#                             outFile.write(chr + "\t")
#                             outFile.write(chrPos + "\t")
#                             outFile.write(str(genomePos) + "\t")
#                             outFile.write("1\t")
#                             outFile.write("\n")
#                             if int(genotypeInfo[2]) >= 10:
#                                 covHomAltCounter = covHomAltCounter + 1
#                             #outFile.write(genotypeInfo[1]+"\t"+genotypeInfo[2])
#         #outFile.write(SNPinfo + "\n")
#         SNParray.append(SNPinfo)
#
# #print("hom="+str(covHomAltCounter)+"\nhet="+str(covHetCounter)+"\n")
# outFile.close()
#
# genomeLen = culumative
#
# genotypeOutName = outName + "_count.txt"
# genotypeOut = open(genotypeOutName, 'w')
# genotypeOut.write("\t\tHeterozygous\tHomozygous\tTotal\t%SNPs_Het\t%Het\tSNPs_per_bp\tHet_SNPs_per_bp\n")
# #genotypeOut.write(str(homozygouteRefCounter)+"\t")
# genotypeOut.write("All_SNPs\t")
# genotypeOut.write(str(totalHetCounter)+"\t")
# genotypeOut.write(str(totalHomAltCounter)+"\t")
# total = totalHomAltCounter+homozygouteRefCounter+totalHetCounter
# genotypeOut.write(str(total)+"\t")
# genotypeOut.write(str((float(totalHetCounter)/total)*100)+"%\t")
# genotypeOut.write(str((float(totalHetCounter)/culumative)*100)+"%\t")
# genotypeOut.write(str('%.3e'%(float(totalHomAltCounter)/genomeLen))+"\t")
# genotypeOut.write(str('%.3e'%(float(totalHetCounter)/genomeLen))+"\n")
# genotypeOut.write("Qual_SNPs\t")
# genotypeOut.write(str(qualHetCounter)+"\t")
# genotypeOut.write(str(qualHomAltCounter)+"\t")
# qualTotal = qualHomAltCounter+homozygouteRefCounter+qualHetCounter
# genotypeOut.write(str(qualHetCounter+qualHomAltCounter)+"\t")
# genotypeOut.write(str((float(qualHetCounter)/qualTotal)*100)+"%\t")
# genotypeOut.write(str((float(qualHetCounter)/culumative)*100)+"%\t")
# genotypeOut.write(str('%.3e'%(float(qualHomAltCounter)/genomeLen))+"\t")
# genotypeOut.write(str('%.3e'%(float(qualHetCounter)/genomeLen))+"\n")
# genotypeOut.write("Cov_SNPs\t")
# genotypeOut.write(str(covHetCounter)+"\t")
# genotypeOut.write(str(covHomAltCounter)+"\t")
# covTotal = covHetCounter+covHomAltCounter
# genotypeOut.write(str(covTotal)+"\t")
# genotypeOut.write(str((float(covHetCounter)/covTotal)*100)+"%\t")
# genotypeOut.write(str((float(covHetCounter)/culumative)*100)+"%\t")
# genotypeOut.write(str('%.3e'%(float(covHomAltCounter)/genomeLen))+"\t")
# genotypeOut.write(str('%.3e'%(float(covHetCounter)/genomeLen))+"\n")
#
