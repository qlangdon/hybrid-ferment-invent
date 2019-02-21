__author__ = 'Quinn'

import Bio
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
import sys, re, time, argparse

parser = argparse.ArgumentParser(description="Call SNPs from whole genome combined FASTA")
parser.add_argument('--out', help="Output prefix, required", required=True)
parser.add_argument('--fasta', help="fasta to use, required", required=True)
parser.add_argument('--chrLength', help="Chromosome length files, required", required=True)
parser.add_argument('--strainKey', help="Key to strains wanted and how to rename them")
parser.add_argument('--N', help="Ignore N threshold, optional")
args = parser.parse_args()

#print(len(sys.argv))
fastaName = args.fasta
prefix = args.out
chrLengthFile = args.chrLength
SNPoutputFile = prefix+"_SNPs.fasta"
SNPposOutputFile = prefix+"_SNPpos.txt"
summaryOutputFile = prefix+"_strainSummary.txt"
ignoreN = float(0.0)
if args.N:
    ignoreN = float(args.N)
else:
    ignoreN = float(0.0)

pickStrains = False
if args.strainKey:
    pickStrains = True
    wantedStrains = list()
    renameDict = {}
    strainKeyFile = open(args.strainKey)
    strainsRename = strainKeyFile.readlines()
    for line in strainsRename:
        currentLine = line.strip('\n')
        strainInfo = currentLine.split('\t')
        strainBase = strainInfo[0]
        newName = strainInfo[1]
        wantedStrains.append(strainBase)
        renameDict[strainBase] = newName


chrLengthDict = {}
lengthFile = open(chrLengthFile, 'r')
chrLengths = lengthFile.read().splitlines()
chrList = ("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")
if len(chrLengths) == 17:
    chrList = ("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrmt")
elif len(chrLengths) == 18:
    chrList = ("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chr2_micron", "chrMT")
counter = 0
for chrLength in chrLengths:
    chrLengthDict[chrList[counter]] = chrLength
    counter += 1

snpDict = {}
strainDict = {}
strainSummary = {}
totalStrainList = list()
seqLength = int()
fasta = open(fastaName, 'r')
for seq_record in SeqIO.parse(fasta, "fasta"):
    strainName = seq_record.id
    totalStrainList.append(strainName)
    strainDict[strainName] = seq_record.seq
    seqLength = len(seq_record.seq)
    snpDict[strainName] = str()
    strainSummary[strainName] = {}
    strainSummary[strainName]['A'] = 0
    strainSummary[strainName]['T'] = 0
    strainSummary[strainName]['C'] = 0
    strainSummary[strainName]['G'] = 0
    strainSummary[strainName]['N'] = 0

if pickStrains is False:
    strainList = totalStrainList
elif pickStrains is True:
    strainList = wantedStrains


SNPpos = str()
first = strainList[0]
refMissing = 0
hasData = 0

SNPposOutput = open(SNPposOutputFile, 'w')
start = time.time()
numSNPs = 0
for i in xrange(0, seqLength): #go through the genome position by position
    compare = strainDict[first][i]
    if re.match('A|T|C|G', compare):
        hasData += 1
        matches = 0
        missing = 0
        SNP = 0
        propMiss = 0.0
        for strain in strainList[1: len(strainList)]:
            if re.match(compare, strainDict[strain][i]):
                matches += 1
            else:
                if re.match('A|T|C|G', strainDict[strain][i]):
                    SNP += 1
                else:
                    missing += 1
        propMiss = float(missing)/len(strainList)
        if SNP>0 and propMiss<=ignoreN:
            numSNPs += 1
            #posOutFile.write("\n"+str(i+1))
            endTime = time.time()-start
            if 60 < endTime < 3600:
                min = int(endTime)/60
                sec = int(endTime-(min*60))
                elapsedTime = str(min) + " mins " + str(sec) + " secs"
            elif 3600 < endTime < 86400:
                hr = int(endTime)/3600
                min = int((endTime - (hr*3600))/60)
                sec = int(endTime - ((hr*60)*60 + (min*60)))
                elapsedTime = str(hr) + " hrs " + str(min) + " mins " + str(sec) + " secs"
            elif 86400 < endTime < 604800:
                day = int(endTime)/86400
                hr = int((endTime-(day*86400))/3600)
                min = int((endTime - (hr*3600+day*86400))/60)
                sec = int(endTime - ((day*86400) + (hr*3600) + (min*60)))
                elapsedTime = str(day)  + " days " + str(hr) + " hrs " + str(min) + " mins " + str(sec) + " secs"
            else:
                elapsedTime = str(int(endTime)) + " secs"
            print("%.5f" % ((i/11633661.0)*100) + "%" + "\tElapsed time: " + elapsedTime + "\tNum SNPs: " + str(numSNPs))
            SNPid = str()
            truPos = i+1
            prevChrEnd = int(0)
            for chr in chrList:
                if truPos >= prevChrEnd:
                    #print(truPos)
                    #print(prevChrEnd)
                    chrPos = truPos-prevChrEnd
                    SNPid = chr + ":" + str(chrPos)
                    prevChrEnd = int(chrLengthDict[chr])
            for strain in strainList:
                nt = strainDict[strain][i]
                if re.match('Y|R|K|M|S|W', nt):
                    nt = "N"
                snpDict[strain] = snpDict[strain] + nt
                strainSummary[strain][nt] += 1
            SNPposOutput.write(str(i)+"\t"+SNPid+"\n")

SNPposOutput.close()

SNPoutput = open(SNPoutputFile, 'w')
for strain in strainList:
    if pickStrains is True:
        newName = renameDict[strain]
        SNPoutput.write(">" + newName + "\n" + snpDict[strain] + "\n")
    else:
        SNPoutput.write(">"+strain+"\n"+snpDict[strain]+"\n")

SNPoutput.close()

summaryOutput = open(summaryOutputFile, 'w')
if pickStrains is True:
    summaryOutput.write("Total SNPs = " + str(numSNPs) + "\nstrainID\trenamedStrain\tperN\tnumN\tperA\tnumA\tperT\tnumT\tperC\tnumC\tperG\tnumG\n")
else :
    summaryOutput.write("Total SNPs = " + str(numSNPs) + "\nstrainID\tperN\tnumN\tperA\tnumA\tperT\tnumT\tperC\tnumC\tperG\tnumG\n")
for strain in strainList:
    N = strainSummary[strain]['N']
    perN = "%.2f" % ((float(N)/numSNPs)*100) + "%"
    A = strainSummary[strain]['A']
    perA = "%.2f" % ((float(A) / numSNPs) * 100) + "%"
    T = strainSummary[strain]['T']
    perT = "%.2f" % ((float(T) / numSNPs) * 100) + "%"
    C = strainSummary[strain]['C']
    perC = "%.2f" % ((float(C) / numSNPs) * 100) + "%"
    G = strainSummary[strain]['G']
    perG = "%.2f" % ((float(G) / numSNPs) * 100) + "%"
    if pickStrains is True:
        newName = renameDict[strain]
        summaryOutput.write(strain + "\t" + newName + "\t"  + perN + "\t" + str(N) + "\t" + perA + "\t" + str(A) + "\t" + perT + "\t" + str(T) + "\t" + perC + "\t" + str(C) + "\t" + perG + "\t" + str(G) + "\n")

    else:
        summaryOutput.write(strain + "\t" + perN + "\t" + str(N) + "\t" + perA + "\t" + str(A) + "\t" + perT + "\t" + str(T) + "\t" + perC + "\t" + str(C) + "\t" + perG + "\t" + str(G) + "\n")

summaryOutput.close()
