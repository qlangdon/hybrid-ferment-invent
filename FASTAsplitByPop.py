__author__ = 'Quinn'

from Bio import SeqIO
import sys, argparse, re, time

parser = argparse.ArgumentParser(description="Take popd strain synonyms to create pop key to fasta name")
parser.add_argument('--fasta', help="Name of fasta to change", required=True)
parser.add_argument('--key', help="pop file", required=True)

args = parser.parse_args()

fastaFile = args.fasta
popFile = args.key

popDict = {}
popList = open(popFile, 'r')
popLines = popList.read().splitlines()
for line in popLines:
    popSplit = line.strip('\n').split('\t')
    pop = popSplit[1]
    strainID = popSplit[0]
    if pop in popDict:
        popDict[pop].append(strainID)
    else:
        popDict[pop] = [strainID]

popList.close()

genomeDict= {}
genomeFASTA = open(fastaFile, 'r')
for seq_record in SeqIO.parse(genomeFASTA, "fasta"):
    strainName = str(seq_record.id)
    strainSeq = str(seq_record.seq)
    genomeDict[strainName] = strainSeq

genomeFASTA.close()

for pop in popDict.keys():
    popStrains = popDict[pop]
    popFASTA = open(pop+".fasta", 'w')
    for strainName in popStrains:
        popFASTA.write(">"+strainName+"\n"+genomeDict[strainName]+"\n")
    popFASTA.close()
