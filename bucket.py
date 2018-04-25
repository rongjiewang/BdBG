from __future__ import print_function

import sys
import screed
import os
import string
import numpy as np
from ngs_plumbing import dna ## DNA to 2bit
from bitstring import * ##glombe integer

from khmer import __version__
from khmer.khmer_args import (sanitize_help, ComboFormatter, info,
                              _VersionStdErrAction)
from readFQ import get_chunks,get_chunksData
from screed import Record
from utils import to_str
import copy

class sortSequenceClass(object):
    """docstring for sortSequence"""
    def __init__(self):
        self.mutiDict = {} ### kmers for buckets 
        self.bucketTable = {} ### buckets index and it's reads number
        self.sequenceTable = {} ####bucket as diction index, within bucket include (seq,reverseFlag,indexPos) 
        self.sequenceTableSave = {}
        self.kmerLen = 8
        self.bucketIndexLen = 8
        self.recordNum = 0
        self.records = [] #read the input data
        self.onsies = 0 ####the number of bucket which only have one seq
        self.leftOnesies = 0 ## the number onesies after re-assigned
        self.sequenceN = []  ### the info of N in read  present by (readID,pos,len)
        self.reversedNum = 0
        self.skipZone = 0
        self.record = Record()
        self.buffer = None
        self.sequenceNIndex = [] #indicate if the seq include N
        self.filePath = ""#indecate the file path
        self.readLen = 0
    def emptyPara(self):
        self.mutiDict = {} ### kmers for buckets 
        self.bucketTable = {} ### buckets index and it's reads number
        self.sequenceTable = {} ####bucket as diction index, within bucket include (seq,reverseFlag,indexPos) 
        self.sequenceTableSave = {}
        self.recordNum = 0
        self.records = [] #read the input data
        self.onsies = 0 ####the number of bucket which only have one seq
        self.leftOnesies = 0 ## the number onesies after re-assigned
        self.sequenceN = []  ### the info of N in read  present by (readID,pos,len)
        self.reversedNum = 0
        self.record = Record()
        self.buffer = None
        self.sequenceNIndex = [] #indicate if the seq include N
    def setskipZone(self,skip):
        self.skipZone = skip
        return
    def getReadLen(self,filename):
        input_iter = screed.open(filename)
        for record in input_iter:
            self.readLen = len(record.sequence)
            break
        return self.readLen
    def setBucketIndex(self,bucketIndexLen):
        self.bucketIndexLen = bucketIndexLen
        return
    def getBucketKmers(self,seq):
        for i in xrange(self.readLen-self.skipZone-self.bucketIndexLen+1):
            yield seq[i:i+self.bucketIndexLen]
    def getReadKmers(self,seq):
        for i in xrange(self.readLen-self.kmerLen+1):
            yield seq[i:i+self.kmerLen]
        
    def analyze_file(self,filename):
        """Run over the given file and count base pairs and sequences."""
        bps = 0
        seqs = 0
        input_iter = screed.open(filename)
        for record in input_iter:
            bps += len(record.sequence)
            seqs += 1
        return bps, seqs
    def getCommenOverlapKmer(self,index,read):
        overlapNum = 0
        if index not in self.mutiDict:
            return overlapNum
        else:
            for kmer in self.getReadKmers(read):
                if kmer in self.mutiDict[index]:
                    overlapNum += 1
            return overlapNum
    def getCommenOverlapKmerNum(self,index,read):
        score = 0
        if index not in self.mutiDict:
            return score
        else:
            for kmer in self.getReadKmers(read):
                if kmer in self.mutiDict[index]:
                    score += self.mutiDict[index][kmer]
            return score
    def getReadScore(self,index,read):
        score = 0
        if index not in self.mutiDict:
            return score
        kmerNumList = []
        for kmer in self.getReadKmers(read):
            if kmer in self.mutiDict[index]:
                kmerNumList.append(self.mutiDict[index][kmer])
        mean = np.mean(kmerNumList)
        score = mean*len(kmerNumList)
        return score

    def addTomutiDict(self,index,read):
        if index in self.mutiDict:
            for kmer in self.getReadKmers(read):
                self.mutiDict[index][kmer] = self.mutiDict[index].get(kmer,0) + 1
        else:
            self.mutiDict.setdefault(index, {})
            for kmer in self.getReadKmers(read):
                self.mutiDict[index][kmer] = self.mutiDict[index].get(kmer,0) + 1
        return
    def getMostCommenIndex(self):
        read = self.record['sequence']
        rcread = screed.rc(str(read))
        freq_kmer = 0
        mostIndex = 0
        indexPos = 0
        reverseFlag = False
        kmerIter = 0
        for kmer in self.getBucketKmers(read):
            m = self.getCommenOverlapKmer(kmer,read)
            if m > freq_kmer:
                freq_kmer = m
                mostIndex = kmer
                indexPos = kmerIter
            kmerIter += 1
        kmerIter = 0
        for kmer in self.getBucketKmers(rcread):
            m = self.getCommenOverlapKmer(kmer,rcread)
            if m > freq_kmer:
                freq_kmer = m
                mostIndex = kmer
                reverseFlag = True
                indexPos = kmerIter
            kmerIter += 1
        #add new kmer to mutiDict
        if freq_kmer > 0:
            if(reverseFlag):
                self.addTomutiDict(mostIndex,rcread)
            else:
                self.addTomutiDict(mostIndex,read)
        else:
            minKmer,reverseFlag,indexPos = self.getReadMinKmers()
            if(reverseFlag):
                self.addTomutiDict(minKmer,rcread)
            else:
                self.addTomutiDict(minKmer,read)

        #find a commen index
        if freq_kmer:
            return mostIndex,reverseFlag,indexPos
        else:
            return minKmer,reverseFlag,indexPos

    def twobit_repr(self,ch):
        x = 0 if ch.upper() == 'A'  else  1 if ch.upper() == 'C' else 2 if ch.upper() == 'G' else 3
        return x               
    def getIndexMinKmer(self,minKmer):
        h = 0
        for i in range(len(minKmer)):
            h = h << 2;
            h |= self.twobit_repr(minKmer[i])
        return h
    def getReadMinKmers(self):
        read = self.record['sequence']
        rcread = screed.rc(str(read))
        minKmer = "T"*self.kmerLen
        reverseFlag = False
        indexPos = 0
        kmerIter = 0
        for kmer in self.getBucketKmers(read[0:(self.readLen-self.skipZone)]):
            if kmer < minKmer:
                minKmer = kmer
                indexPos = kmerIter
            kmerIter +=1

        kmerIter = 0
        for kmer in self.getBucketKmers(rcread[0:(self.readLen-self.skipZone)]):
            if kmer < minKmer:
                minKmer = kmer
                reverseFlag = True
                indexPos = kmerIter
            kmerIter += 1
        return minKmer,reverseFlag,indexPos
    def sortSeqence(self):
        """Run over the given file and sort the sequences."""
        index,reverseFlag,indexPos = self.getMostCommenIndex()
        self.bucketTable[index] = self.bucketTable.get(index,0) + 1
        if reverseFlag:
            self.record['sequence'] = screed.rc(str(self.record['sequence']))
        self.sequenceTable.setdefault(index, []).append([copy.deepcopy(self.record),reverseFlag,indexPos])
    def reassigned(self): 
        for countIndex in self.bucketTable.keys():
            if self.bucketTable[countIndex] == 1:
                del self.mutiDict[countIndex]
                del self.bucketTable[countIndex]
                record = self.sequenceTable[countIndex][0]
                del self.sequenceTable[countIndex]
                self.record = record[0]
                index,reverseFlag,indexPos = self.getMostCommenIndex()                  
                self.bucketTable[index] = self.bucketTable.get(index,0) + 1
                self.sequenceTable.setdefault(index, []).append([record[0],reverseFlag,indexPos])
   
    def outputInfo(self):
        print("the bucket index Length is ",self.bucketIndexLen)
        print("the kmer index Length is ",self.kmerLen)
        for countIndex in self.bucketTable.keys():
            if self.bucketTable[countIndex] == 1:
                self.onsies += 1
        print("the number onsies is",self.onsies)
        print("the number of bucket is",len(self.bucketTable.keys()))
        print("the number of read is",self.recordNum)
        print("the number of reversed reads is {0},occupy {1}% ".format(self.reversedNum, self.reversedNum*100/self.recordNum))
    def openInput(self,fileName):
        self.records = screed.open(fileName)

    def sort_key(self,seq):
        return seq[2],seq[0].sequence[(seq[2]+self.kmerLen-1):]

    def setKmerLen(self,L):
        self.kmerLen = L
        return
    def setOutPath(self,path):
        self.filePath = path
        return
    def initialFile(self):
        fileoutBucketIndex    =  self.filePath + ".index"
        fileoutBucketCov    =  self.filePath + ".cov"
        fileoutRCBin        = self.filePath +".rc"
        fileoutSequenceN = self.filePath +".N"
        fileoutIndexPos  = self.filePath +".indexPos"
                            
        if os.path.exists(fileoutBucketIndex):
            os.remove(fileoutBucketIndex)
        if os.path.exists(fileoutBucketCov):
            os.remove(fileoutBucketCov)      
        if os.path.exists(fileoutRCBin):
            os.remove(fileoutRCBin) 
        if os.path.exists(fileoutSequenceN):
            os.remove(fileoutSequenceN)
        if os.path.exists(fileoutIndexPos):
            os.remove(fileoutIndexPos) 
    def outPutSeqence(self):
        fileoutBucketIndex    =  self.filePath + ".index"
        fileoutBucketCov    =  self.filePath + ".cov"
        fileoutRCBin        = self.filePath +".rc"
        fileoutSequenceN = self.filePath +".N"
        fileoutIndexPos  = self.filePath +".indexPos"
                            
        BucketIndexFile = open(fileoutBucketIndex,"a+")
        BucketCovFile = open(fileoutBucketCov,"a+")        
        RCBinFile = open(fileoutRCBin,"a+")       
        SequenceNFile = open(fileoutSequenceN,"a+")
        IndexPosFile = open(fileoutIndexPos,"a+") 

        ###output bucket########
        curren_index = 0
        for index in sorted(self.bucketTable):
            indexBin = self.getIndexMinKmer(index) -curren_index
            curren_index = indexBin + curren_index
            BucketIndexFile.write(str(indexBin)+"\n")###store index in number
            BucketCovFile.write(str(self.bucketTable[index]-1)+"\n")###store covrage in number
        BucketIndexFile.close()
        BucketCovFile.close()

        ###output sequence#######
        IDnum = 1
        RCbin = BitStream()

        for index in sorted(self.sequenceTable.keys()):
            prePos = 0  ###for delta encode index position
            for record in sorted(self.sequenceTable[index],key=self.sort_key):
                if(len(self.sequenceTable[index])==1):
                    next
                if record[1]:
                    self.reversedNum += 1
                IDnum += 1
                RCbin.append('0b1' if record[1] else '0b0')
                IndexPosFile.write(chr(record[2]-prePos))
                prePos = record[2]
        RCbin.tofile(RCBinFile)


        RCBinFile.close()
        IndexPosFile.close()

        ####output N pos and Length information#######
        Nindex = BitStream()
        for Nflag in self.sequenceNIndex:
            Nindex.append('0b1' if Nflag else '0b0')
        Nindex.tofile(SequenceNFile)
        for Ninfo in self.sequenceN:
            prePos = 0
            for pos in Ninfo[0]:
                SequenceNFile.write(chr(pos-prePos))
                prePos = pos
            for posLen in Ninfo[1]:
                SequenceNFile.write(chr(posLen-1))
            SequenceNFile.write("\n")
        SequenceNFile.close()
        return
    def getNInRead(self):
        read = self.record['sequence']
        pos = []
        posLen = []
        i = 0        
        while(i < len(read)):
            if read[i] == "N":
                nLen = 1
                pos.append(i)
                i += 1
                while (i < len(read) ):
                    if read[i] == "N":
                        nLen += 1
                        i += 1
                    else:
                        break
                posLen.append(nLen)
            i += 1
        return (pos,posLen)

    def replaceN(self):
        if "N" not in self.record['sequence']:
            self.sequenceNIndex.append(0)
        else:
            self.sequenceNIndex.append(1)
            pos,length = self.getNInRead()
            self.sequenceN.append([pos,length])
            self.record['sequence'] = string.replace(self.record['sequence'], "N", "T")
        return

    def fastq_iter(self, line=None, parse_description=False):
        """
        Iterator over the given FASTQ file handle returning records. buffer
        is a handle to a data set
        """
        if line is None:
            line = self.buffer.next()
        line = to_str(line.strip())
        while line:
            if line and not line.startswith('@'):
                raise IOError("Bad FASTQ format: no '@' at beginning of line")

            # Try to grab the name and (optional) annotations
            if parse_description:
                try:
                    self.record['name'], self.record['annotations'] = line[1:].split(' ', 1)
                except ValueError:  # No optional annotations
                    self.record['name'] = line[1:]
                    self.record['annotations'] = ''
                    pass
            else:
                self.record['name'] = line[1:]
                self.record['annotations'] = ''

            # Extract the sequence lines
            sequence = []
            line = to_str(self.buffer.next())
            while line and not line.startswith('+') and not line.startswith('#'):
                sequence.append(line)
                line = to_str(self.buffer.next())

            self.record['sequence'] = ''.join(sequence)

            # Extract the quality lines
            quality = []
            line = to_str(self.buffer.next())
            seqlen = len(self.record['sequence'])
            aclen = 0
            while not line == '' and aclen < seqlen:
                quality.append(line)
                aclen += len(line)
                line = to_str(self.buffer.next())

            self.record['quality'] = ''.join(quality)
            if len(self.record['sequence']) != len(self.record['quality']):
                raise IOError('sequence and quality strings must be '
                              'of equal length')
            yield self.record
    def setBuffer(self,buffer):
        self.buffer = iter(buffer)
        return
    def saveSeqTable(self):
        for bucketIndex in sorted(self.sequenceTable.keys()):
            for record in self.sequenceTable[bucketIndex]:
                self.sequenceTableSave.setdefault(bucketIndex, []).append([str(record[0]['sequence']),record[2]])
        return
    def compressFile(self):
        fileoutBucketIndex    =  self.filePath + ".index"
        fileoutBucketCov    =  self.filePath + ".cov"
        fileoutSequenceN = self.filePath +".N"
        fileoutIndexPos  = self.filePath +".indexPos"

        if os.path.exists(fileoutBucketIndex + ".lz"):
            os.remove(fileoutBucketIndex+ ".lz")
        if os.path.exists(fileoutBucketCov+ ".lz"):
            os.remove(fileoutBucketCov+ ".lz")
        if os.path.exists(fileoutSequenceN+ ".lz"):
            os.remove(fileoutSequenceN+ ".lz")
        if os.path.exists(fileoutIndexPos+ ".lz"):
            os.remove(fileoutIndexPos+ ".lz")
        comand= "plzip " + fileoutBucketIndex
        os.system(comand)
        comand= "plzip " + fileoutBucketCov
        os.system(comand)
        comand= "plzip " + fileoutSequenceN
        os.system(comand)
        comand= "plzip " + fileoutIndexPos
        os.system(comand)
    def analysisFileSize(self):
        fileoutBucketIndex    =  self.filePath + ".index.lz"
        fileoutBucketCov    =  self.filePath + ".cov.lz"
        fileoutSequenceN = self.filePath +".N.lz"
        fileoutIndexPos  = self.filePath +".indexPos.lz"
        fileoutRCBin = self.filePath + ".rc"
        indexsize=os.path.getsize(fileoutBucketIndex)
        covsize=os.path.getsize(fileoutBucketCov)
        possize=os.path.getsize(fileoutIndexPos)
        Nsize=os.path.getsize(fileoutSequenceN)
        RCsize=os.path.getsize(fileoutRCBin)
        filesize=[indexsize,covsize,possize,Nsize,RCsize]
        print("[index,cov,pos,N,RC]")
        print(filesize,"sum",sum(filesize),"bytes")




