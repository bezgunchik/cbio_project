import sys
import os

scoreDict = dict()


def unitReadScore(read, start):
    ASCII = read[3][int(start) - int(read[1])]
    phredQScore = ord(ASCII) - 33
    error = 10 ** (-phredQScore/10)
    return error


def getSeqQ(snpName, chr, start, end, bamFile):


    command = "samtools view -F 1796 -q 10 " \
              f"{bamFile} " \
              f"{chr}:{start}-{end} " \
              "| grep '150M' | awk '{print $3,$4,$10,$11}'"
    stream = os.popen(command)
    # print(stream.read())
    readsNum = 0
    currentScore = 0
    for read in stream.readlines():
        currentScore += unitReadScore(read.split(), start)
        readsNum += 1
    scoreDict[snpName] = currentScore / readsNum
    print('Average sequencing error of ' + snpName + ': ' + str(scoreDict[snpName]))
    #samtools view -F 1796 -q 10 CNVS-NORM-110033745-cfDNA-WGBS-Rep1.bam chr3:11620414-11620415 | grep '150M' | head -1 | awk '{print $3,$4,$10,$11}'



if __name__ == '__main__':
    args = sys.argv
    with open(args[1], 'r') as snpFile:
        for line in snpFile:
            splitLine = line.split('\t')
            snpName = splitLine[3][:-1]
            scoreDict[snpName] = 0
            chr = splitLine[0]
            start = splitLine[1]
            end = splitLine[2]
            getSeqQ(snpName, chr, start, end, args[2])














