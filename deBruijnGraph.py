from __future__ import print_function
from __future__ import division

import sys
sys.path.append("./arithCode")
import arithmeticcoding
import os
from ngs_plumbing import dna ## DNA to 2bit
from bitstring import * ##glombe integer
from collections import defaultdict
import dna #for reverse complement dna
import numpy as np
import bitio
from bitstream import BitStream as sream


def multi_dimensions(n, type):
  """ Creates an n-dimension dictionary where the n-th dimension is of type 'type'
  """  
  if n<=1:
    return type()
  return defaultdict(lambda:multi_dimensions(n-1, type))
def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

class encodeGraphClass(object):
    """docstring for sortSequence"""
    def __init__(self,indexLen,kmerLen,path,bucket_dict={},sequenceTable={}):
        self.mutiDict = bucket_dict ### kmers for buckets 
        self.sequenceTable = sequenceTable
        self.kmerLen = kmerLen
        self.indexLen = indexLen 
        #self.bucketDict = defaultdict(lambda : defaultdict(dict))
        self.bucketDict = {}#nested_dict(2, int)
        self.encodeBucketPath = {}
        self.newNodeNum = 0
        self.simpleNodeNum = 0
        self.tipNodeNum = 0
        self.bifurNodeNum = 0
        self.deleteBifurRatio = 0.2
        self.outPutPath = path
        self.dna2num={"A":0,"C":1,"G":2,"T":3}
        self.num2dna={0:"A",1:"C",2:"G",3:"T"}
        self.dna2bit={"A":'0b00',"C":'0b01',"G":'0b10',"T":'0b11'}
        self.num2dna={0:"A",1:"C",2:"G",3:"T"}
        self.firstSeq = BitStream()
        self.numFlag = BitStream()
        self.freq3 =  arithmeticcoding.SimpleFrequencyTable(arithmeticcoding.FlatFrequencyTable(3))
        self.freq4 =  arithmeticcoding.SimpleFrequencyTable(arithmeticcoding.FlatFrequencyTable(4))
        self.freqs =  arithmeticcoding.SimpleFrequencyTable(arithmeticcoding.FlatFrequencyTable(4))
        self.bitoutL = arithmeticcoding.BitOutputStream(open(self.outPutPath+".bifurL", "wb"))
        self.bitoutR = arithmeticcoding.BitOutputStream(open(self.outPutPath+".bifurR", "wb"))
        self.encodeSeqPathL = self.openFileLeft()
        self.encodeSeqPathR = self.openFileRight()
    def openFileLeft(self):
        enc  = arithmeticcoding.ArithmeticEncoder(self.bitoutL)
        return enc

    def openFileRight(self):
        enc = arithmeticcoding.ArithmeticEncoder(self.bitoutR)
        return enc 
    def emptyPara(self):
        self.mutiDict.clear() ### kmers for buckets 
        self.sequenceTable.clear()
        self.bucketDict.clear()#nested_dict(2, int)
        self.encodeBucketPath.clear()
        self.newNodeNum = 0
        self.simpleNodeNum = 0
        self.tipNodeNum = 0
        self.bifurNodeNum = 0
        self.firstSeq = BitStream()
        self.numFlag = BitStream()
        return
    def setDiction(self,dict):
        self.mutiDict = dict
    def setSequenceTable(self,table):
        self.sequenceTable = table
    def setOutPath(self,path):
        self.outPutPath = path
        return
    def initialOutFile(self):
        fileoutFirstSeq = self.outPutPath + ".firSeq"
        fileoutNumFlag = self.outPutPath + ".numFlag"
        if os.path.exists(fileoutFirstSeq):
            os.remove(fileoutFirstSeq)
        if os.path.exists(fileoutNumFlag):
            os.remove(fileoutNumFlag)
    def outputEncodePath(self):
        self.encodeSeqPathL.finish()
        self.encodeSeqPathR.finish()
        self.bitoutL.close()
        self.bitoutR.close()
        fileoutFirstSeq = self.outPutPath + ".firSeq"
        FirstSeqFile = open(fileoutFirstSeq,"a+")
        self.firstSeq.tofile(FirstSeqFile)
        fileoutNumFlag = self.outPutPath + ".numFlag"
        numFlagFile = open(fileoutNumFlag,"a+")
        self.numFlag.tofile(numFlagFile)
    def twobit_repr(self,ch):
        x = 0 if ch.upper() == 'A'  else  1 if ch.upper() == 'C' else 2 if ch.upper() == 'G' else 3
        return x
    def isDNA(self,ch):
        if ch.upper() == 'A'  or ch.upper() == 'C' or ch.upper() == 'G' or ch.upper() == 'T':
            return True
        return False        
    def analyKmerScore(self):
        for bucketIndex in self.sequenceTable.keys():
            i = 0
            for seq,pos in self.sequenceTable[bucketIndex]:
                score = self.getSeqScore(bucketIndex,seq)
                self.sequenceTable[bucketIndex][i].append(score)
                i += 1
        return
    def sortSeqInbucket(self):
        for bucketIndex in self.sequenceTable.keys():
            self.sequenceTable[bucketIndex] = sorted(self.sequenceTable[bucketIndex],key=self.sort_key, reverse=True)
    def sort_key(self,seq):
        return seq[1],seq[2]
    def getSeqScore(self,bucketIndex,seq):
        kmerNumList = []
        for kmer in self.getKmer(seq):
            if kmer in self.mutiDict[bucketIndex]:
                kmerNumList.append(self.mutiDict[bucketIndex][kmer])
        score = np.sum(kmerNumList)
        return score

    def getKmer(self,seq):
        for i in xrange(len(seq)-self.kmerLen+1):
            yield seq[i:i+self.kmerLen]
    def getKmerR(self,seq,pos):
        for i in xrange(pos,len(seq)-self.kmerLen):
            yield seq[i:i+self.kmerLen+1]
    def getKmerL(self,seq,pos):
        for i in xrange(pos,0,-1):
            yield seq[i-1:i+self.kmerLen]
    def encodeBucket(self):
        seqNum = 1
        for bucketIndex in sorted(self.sequenceTable.keys()):
            iterSeq = iter(self.sequenceTable[bucketIndex])
            seq, pos = next(iterSeq)
            self.encodeFirSeq(seq, pos)
            self.addtoGraph(seq, pos)
            for seq,pos in iterSeq:
                self.encodeSeq(bucketIndex, seq, pos)
                self.addtoGraph(seq,pos)
                if seqNum % 10000 == 0:
                    print('process %d records'% seqNum)
                seqNum += 1
            del self.sequenceTable[bucketIndex]
            self.bucketDict.clear()
        return
    def isUniqueSuffix(self,bucketIndex,kmer):
        for alterKmer in self.getAlterKmerInLastPos(kmer):
            if alterKmer in self.bucketDict.keys():
                return False
        return True
    def isUniquePrefix(self,bucketIndex,kmer):
        for alterKmer in self.getAlterKmerInFirstPos(kmer):
            if alterKmer in self.bucketDict.keys():
                return False
        return True
    def getAlterKmerInLastPos(self,kmer):
        for base in "ACGT":
            if base != kmer[-1]:
                yield kmer[:-1]+base
    def getAlterKmerInFirstPos(self,kmer):
        for base in "ACGT":
            if base != kmer[0]:
                yield base + kmer[1:]
    def Suffix(self,kmer):
        return kmer[1:]
    def Successors_in_graph(self, K):
        succ = []
        if K in self.bucketDict:
            for base in range(4):
                if self.bucketDict[K][base] > 1 and self.sufficient(K,base):
                    succ.append(K[1:]+self.num2dna[base])
        return succ
    def sufficient(self,K,base):
        totalSum = sum(self.bucketDict[K])
        if self.bucketDict[K][base]/totalSum < self.deleteBifurRatio:
            return False
        return True
    def encodeFirSeq(self,seq,pos):# encode the first sequence in bucket
        for base in seq[:pos]:
            self.firstSeq.append(self.dna2bit[base])
        for base in seq[pos+self.kmerLen:]:
            self.firstSeq.append(self.dna2bit[base])
        return
    def encodeSeq(self,bucketIndex,seq,pos):
        #encode the left part
        K = dna.reverse_complement(str(bucketIndex))
        for base in dna.reverse_complement(seq[:pos]):
            pred = self.Successors_in_graph(K)
            if len(pred) == 1:
                if pred[0][-1] != base:
                    self.newNodeNum += 1
                    self.numFlag.append('0b1')
                    symbol = self.dna2num[base]
                    if base > pred[0][-1]:
                        symbol = self.dna2num[base] - 1
                    self.encodeSeqPathL.write(self.freq3,symbol)#save the reverse complement sequence
                else:
                    self.simpleNodeNum += 1
                    self.numFlag.append('0b0')
                K = pred[0]
            else:
                if len(pred) == 0:
                    self.tipNodeNum += 1
                    self.encodeSeqPathL.write(self.freq4,self.dna2num[base])
                else:
                    self.bifurNodeNum += 1
                    self.getFreqs(K)
                    self.encodeSeqPathL.write(self.freqs,self.dna2num[base])
                K = self.Suffix(K) + base
        #encode the right part
        K = str(bucketIndex)
        for base in seq[pos+self.indexLen:]:
            succ = self.Successors_in_graph(K)
            if len(succ) == 1:
                if succ[0][-1] != base:
                    self.newNodeNum += 1
                    self.numFlag.append('0b1')
                    symbol = self.dna2num[base]
                    if base > succ[0][-1]:
                        symbol = self.dna2num[base] - 1
                    self.encodeSeqPathR.write(self.freq3, symbol)
                else:
                    self.simpleNodeNum += 1
                    self.numFlag.append('0b0')
                K = succ[0]
            else:
                if len(succ) == 0:
                    self.tipNodeNum += 1
                    self.encodeSeqPathR.write(self.freq4, self.dna2num[base])
                else:
                    self.bifurNodeNum += 1
                    self.getFreqs(K)
                    self.encodeSeqPathR.write(self.freqs, self.dna2num[base])
                K = self.Suffix(K) + base
        return
    def getFreqs(self,K):
        for base in range(4):
            self.freqs.set(base,int(self.bucketDict[K][base]))
    def outputNodeInfo(self):
        if (self.simpleNodeNum + self.newNodeNum + self.tipNodeNum + self.bifurNodeNum) == 0:
            print("All the node is Tip node")
        else:
            print("the numbers of simpleNodeNum is %d, take up %.2f%%"%(self.simpleNodeNum,self.simpleNodeNum*100.0/(self.simpleNodeNum+\
                  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
            print("the numbers of newNodeNum is %d, take up %.2f%%"%(self.newNodeNum,self.newNodeNum*100.0/(self.simpleNodeNum+\
                  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
            print("the numbers of tipNodeNum is %d, take up %.2f%%"%(self.tipNodeNum,self.tipNodeNum*100.0/(self.simpleNodeNum+\
                  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
            print("the numbers of bifurNodeNum is %d, take up %.2f%%"%(self.bifurNodeNum,self.bifurNodeNum*100.0/(self.simpleNodeNum+\
                  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
    def addtoGraph(self,seq,pos):
        for kmer in self.getKmerR(seq,pos):
            if str(kmer[:-1]) not in self.bucketDict:
              self.bucketDict.setdefault(str(kmer[:-1]), [1,1,1,1])
            self.bucketDict[str(kmer[:-1])][self.dna2num[str(kmer[-1])]] += 1
        for kmer in self.getKmerL(seq,pos):
            if dna.reverse_complement(kmer[1:]) not in self.bucketDict:
                self.bucketDict.setdefault(dna.reverse_complement(kmer[1:]), [1,1,1,1])           
            self.bucketDict[dna.reverse_complement(kmer[1:])][3-self.dna2num[kmer[0]]] += 1
        return
    def outputBucketDiction(self):
        fileoutBucketDictionPath = self.outPutPath + "bucketDiction"
        if os.path.exists(fileoutBucketDictionPath):
            os.remove(fileoutBucketDictionPath)
        BucketDictionPathFile = open(fileoutBucketDictionPath,"w+")
        for bucketIndex in sorted(self.bucketDict.keys()):
            BucketDictionPathFile.write("bucketIndex is "+ bucketIndex+'\n')
            for kmerIndex in self.bucketDict.keys():
                BucketDictionPathFile.write(kmerIndex+'\t'+str(self.bucketDict[kmerIndex])\
                                            +'\n')
        BucketDictionPathFile.close()
        return
    def compressFile(self):
        if os.path.exists(self.outPutPath + ".firSeq.lz"):
            os.remove(self.outPutPath + ".firSeq.lz")
        os.system("plzip " + self.outPutPath + ".firSeq")
        if os.path.exists(self.outPutPath + ".numFlag.lz"):
            os.remove(self.outPutPath + ".numFlag.lz")
        os.system("plzip " + self.outPutPath + ".numFlag")
        return
    def analysisFileSize(self):
        leftsize=os.path.getsize(self.outPutPath + ".bifurL")
        rightsize=os.path.getsize(self.outPutPath + ".bifurR")
        firstSeqsize = os.path.getsize(self.outPutPath + ".firSeq.lz")
        numFlagsize = os.path.getsize(self.outPutPath + ".numFlag.lz")
        filesize=[leftsize,rightsize,firstSeqsize,numFlagsize]
        print("[bifurL,bifurR,firSeq,numFlag]")
        print(filesize,"sum",sum(filesize),"bytes")
        return

class decodeGraphClass(object):
    """docstring for sortSequence"""
    def __init__(self,indexLen,kmerLen,inPutPath, outPutPath, bucket_dict={},sequenceTable={}):
        self.mutiDict = bucket_dict ### kmers for buckets 
        self.sequenceTable = sequenceTable
        self.kmerLen = kmerLen
        self.indexLen = indexLen 
        self.bucketDict = {}#nested_dict(2, int)
        self.encodeBucketPath = {}
        self.newNodeNum = 0
        self.simpleNodeNum = 0
        self.tipNodeNum = 0
        self.bifurNodeNum = 0
        self.deleteBifurRatio = 0.2
        self.inPutPath = inPutPath
        self.outPutPath = outPutPath
        self.dna2num={"A":0,"C":1,"G":2,"T":3}
        self.num2dna={0:"A",1:"C",2:"G",3:"T"}
        self.recdna={"A":"T","C":"G","G":"C","T":"A"}        
        self.dna2bit={"A":'0b00',"C":'0b01',"G":'0b10',"T":'0b11'}
        self.num2dna={0:"A",1:"C",2:"G",3:"T"}
        self.firstSeq = BitStream()
        self.numFlag = BitStream()
        self.freq3 =  arithmeticcoding.SimpleFrequencyTable(arithmeticcoding.FlatFrequencyTable(3))
        self.freq4 =  arithmeticcoding.SimpleFrequencyTable(arithmeticcoding.FlatFrequencyTable(4))
        self.freqs =  arithmeticcoding.SimpleFrequencyTable(arithmeticcoding.FlatFrequencyTable(4))
        self.bitoutL = arithmeticcoding.BitInputStream(open(self.inPutPath+".bifurL", "rb"))
        self.bitoutR = arithmeticcoding.BitInputStream(open(self.inPutPath+".bifurR", "rb"))
        self.decodeSeqPathL = self.openFileLeft()
        self.decodeSeqPathR = self.openFileRight()
        self.outFileName = outPutPath + ".dna"
        self.outFile = self.openOutFile()
        self.readLen = 100
        self.seqence = ""
        self.bucketIndex = [] #bucket index 
        self.bucketCov = [] # reads number in bucket
        self.readIndexPos = [] #index positions in each read 
        self.readrc = sream() # read in forward or backward
        self.readN = {"flag":sream(), "pos":[], "l":[]} # N in read indicate, number, position and length
        self.numFlag = sream() #new nodes indicate
    def setMetaInfor(self, readLen, bucketIndexLen):
        self.readLen = readLen
        self.indexLen = bucketIndexLen
        self.kmerLen = bucketIndexLen
        return
    def loadBucktData(self,bucketIndex, bucketCov, readIndexPos, readrc, readN):
        self.bucketIndex = bucketIndex #bucket index 
        self.bucketCov = bucketCov # reads number in bucket
        self.readIndexPos = readIndexPos #index positions in each read
        self.readrc = readrc # read in forward or backward
        self.readN = readN # N in read indicate, position and length
    def openFileLeft(self):
        enc  = arithmeticcoding.ArithmeticDecoder(self.bitoutL)
        return enc
    def openOutFile(self):
        self.initialOutFile()
        enc = open(self.outFileName,'a+')
        return enc

    def openFileRight(self):
        enc = arithmeticcoding.ArithmeticDecoder(self.bitoutR)
        return enc 
    def setDiction(self,dict):
        self.mutiDict = dict
    def setSequenceTable(self,table):
        self.sequenceTable = table
    def setPath(self,inPutPath, outPutPath):
        self.inPutPath = inPutPath
        self.outPutPath = outPutPath
        return
    def initialOutFile(self):
        if os.path.exists(self.outFileName):
            os.remove(self.outFileName)
    def twobit_repr(self,ch):
        x = 0 if ch.upper() == 'A'  else  1 if ch.upper() == 'C' else 2 if ch.upper() == 'G' else 3
        return x
    def isDNA(self,ch):
        if ch.upper() == 'A'  or ch.upper() == 'C' or ch.upper() == 'G' or ch.upper() == 'T':
            return True
        return False        
    def analyKmerScore(self):
        for bucketIndex in self.sequenceTable.keys():
            i = 0
            for seq,pos in self.sequenceTable[bucketIndex]:
                score = self.getSeqScore(bucketIndex,seq)
                self.sequenceTable[bucketIndex][i].append(score)
                i += 1
        return
    def sortSeqInbucket(self):
        for bucketIndex in self.sequenceTable.keys():
            self.sequenceTable[bucketIndex] = sorted(self.sequenceTable[bucketIndex],key=self.sort_key, reverse=True)
    def sort_key(self,seq):
        return seq[1],seq[2]
    def getSeqScore(self,bucketIndex,seq):
        kmerNumList = []
        for kmer in self.getKmer(seq):
            if kmer in self.mutiDict[bucketIndex]:
                kmerNumList.append(self.mutiDict[bucketIndex][kmer])
        score = np.sum(kmerNumList)
        return score

    def getKmer(self,seq):
        for i in xrange(len(seq)-self.kmerLen+1):
            yield seq[i:i+self.kmerLen]
    def getKmerR(self,seq,pos):
        for i in xrange(pos,len(seq)-self.kmerLen):
            yield seq[i:i+self.kmerLen+1]
    def getKmerL(self,seq,pos):
        for i in xrange(pos,0,-1):
            yield seq[i-1:i+self.kmerLen]
    def isUniqueSuffix(self,bucketIndex,kmer):
        for alterKmer in self.getAlterKmerInLastPos(kmer):
            if alterKmer in self.bucketDict.keys():
                return False
        return True
    def isUniquePrefix(self,bucketIndex,kmer):
        for alterKmer in self.getAlterKmerInFirstPos(kmer):
            if alterKmer in self.bucketDict.keys():
                return False
        return True
    def getAlterKmerInLastPos(self,kmer):
        for base in "ACGT":
            if base != kmer[-1]:
                yield kmer[:-1]+base
    def getAlterKmerInFirstPos(self,kmer):
        for base in "ACGT":
            if base != kmer[0]:
                yield base + kmer[1:]
    def Suffix(self,kmer):
        return kmer[1:]
    def Successors_in_graph(self, K):
        succ = []
        if K in self.bucketDict:
            for base in range(4):
                if self.bucketDict[K][base] > 1 and self.sufficient(K,base):
                    succ.append(K[1:]+self.num2dna[base])
        return succ
    def sufficient(self,K,base):
        totalSum = sum(self.bucketDict[K])
        if self.bucketDict[K][base]/totalSum < self.deleteBifurRatio:
            return False
        return True
    def getFreqs(self,K):
        for base in range(4):
            self.freqs.set(base,int(self.bucketDict[K][base]))
    def addtoGraph(self,seq,pos):
        for kmer in self.getKmerR(seq,pos):
            if str(kmer[:-1]) not in self.bucketDict:
              self.bucketDict.setdefault(str(kmer[:-1]), [1,1,1,1])
            self.bucketDict[str(kmer[:-1])][self.dna2num[str(kmer[-1])]] += 1
        for kmer in self.getKmerL(seq,pos):
            if dna.reverse_complement(kmer[1:]) not in self.bucketDict:
                self.bucketDict.setdefault(dna.reverse_complement(kmer[1:]), [1,1,1,1])           
            self.bucketDict[dna.reverse_complement(kmer[1:])][3-self.dna2num[kmer[0]]] += 1
        return
    def decompress(self):
        self.unPlzipFiles()
        self.decodeSeq()
        self.outFile.close()
        return
    def unPlzipFiles(self):
        for name in ["firSeq","numFlag"]:
            file = self.inPutPath + "." + name + ".lz"                           
            if os.path.exists(file):
                comand= "plzip -d " + file
                os.system(comand)
        return
    def decodeSeq(self):
        self.decodeNumFlag()
        file = self.inPutPath + ".firSeq"
        firSeqFile = bitio.bit_open(file, "r")
        firSeqLoadLen = (self.readLen-self.indexLen)*2
        for bucketIndex in self.bucketIndex:
            self.bucketDict.clear() #initial the bucketDiction with nothing
            indexPos = self.readIndexPos.pop(0)
            self.sequence = self.decompressFirSeq(firSeqFile, bucketIndex, indexPos, firSeqLoadLen)
            self.addtoGraph(self.sequence, indexPos)
            self.outPutSeqence()
            readNum = 0
            bucketReads = self.bucketCov.pop(0)
            while readNum < bucketReads:
                indexPos += self.readIndexPos.pop(0)
                self.sequence = self.decompressGraph(bucketIndex, indexPos)
                self.addtoGraph(self.sequence, indexPos)
                self.outPutSeqence()
                readNum += 1
        print("simpleNodeNum is %d"% self.simpleNodeNum)
        print("newNodeNum is %d"% self.newNodeNum)
        print("tipNodeNum is %d"% self.tipNodeNum)
        print("bifurNodeNum is %d"% self.bifurNodeNum)
        firSeqFile.close()
        self.outFile.close()
        print("decoding success!")
        return
    def decompressFirSeq(self, firSeqFile, bucketIndex, indexPos, firSeqLoadLen):
        seq = firSeqFile.read_bits(firSeqLoadLen)
        sequence = self.num2Base(seq)
        sequence.insert(indexPos, bucketIndex)
        sequence = ''.join(sequence)
        return sequence
    def num2Base(self, indexNum):
        seq = []
        while indexNum:
            seq.append(self.num2dna[(indexNum & 3)])
            indexNum = indexNum >> 2
        while len(seq) < (self.readLen-self.indexLen):
            seq.append("A")
        return seq[::-1]
    def outPutSeqence(self):
        if self.readrc.read(bool,1)[0]:
            self.sequence = dna.reverse_complement(self.sequence)
        self.replaceN()        
        self.outFile.write(self.sequence+"\n")
        return
    def replaceN(self):
        if self.readN["flag"].read(bool,1)[0]:
            prePos = 0
            for pos, length in zip(self.readN["pos"].pop(0), self.readN["l"].pop(0)):
                pos += prePos
                self.sequence = self.sequence[:pos] + "N"*int(length+1) + self.sequence[(pos+length+1):]
                prePos = pos
        return
    def decompressGraph(self,bucketIndex, indexPos):
        sequence = bucketIndex
        #decode the left path
        K = dna.reverse_complement(bucketIndex)
        base = ""
        for i in range(indexPos):
            pred = self.Successors_in_graph(K)
            if len(pred) == 1:
                if self.numFlag.read(bool,1)[0]: # New node:
                    self.newNodeNum += 1
                    symbol = self.decodeSeqPathL.read(self.freq3)
                    base = self.num2dna[symbol]
                    if base >= pred[0][-1]:
                        base = self.num2dna[symbol+1]
                else:
                    self.simpleNodeNum += 1
                    base = pred[0][-1]
                K = pred[0]
            else:
                if len(pred) == 0:
                    self.tipNodeNum += 1
                    symbol = self.decodeSeqPathL.read(self.freq4)
                    base = self.num2dna[symbol]
                else:
                    self.bifurNodeNum += 1
                    self.getFreqs(K)
                    symbol = self.decodeSeqPathL.read(self.freqs)
                    base = self.num2dna[symbol]
                K = self.Suffix(K) + base
            sequence = self.recdna[base] + sequence
        #decode the right path
        K = bucketIndex
        for i in range(self.readLen-indexPos-self.indexLen):
            succ = self.Successors_in_graph(K)
            if len(succ) == 1:
                if self.numFlag.read(bool,1)[0]: # New node
                    self.newNodeNum += 1
                    symbol = self.decodeSeqPathR.read(self.freq3)
                    base = self.num2dna[symbol]
                    if base >= succ[0][-1]:
                        base = self.num2dna[symbol+1]  
                else:
                    self.simpleNodeNum += 1
                    base = succ[0][-1]
                K = succ[0]
            else:
                if len(succ) == 0:
                    self.tipNodeNum += 1
                    symbol = self.decodeSeqPathR.read(self.freq4)
                    base = self.num2dna[symbol]
                else:
                    self.bifurNodeNum += 1
                    self.getFreqs(K)
                    symbol = self.decodeSeqPathR.read(self.freqs)
                    base = self.num2dna[symbol]
                K = self.Suffix(K) + base
            sequence += base
        return sequence
    def decodeNumFlag(self):
        file = self.inPutPath + ".numFlag"
        with open(file,'rb') as f:
            data=f.read()
            self.numFlag.write(data,bytes)
        f.close()
        return




