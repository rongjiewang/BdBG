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

class graphClass(object):
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
            firSeqFlag = True
            for seq,pos in self.sequenceTable[bucketIndex]:
                self.encodeSeq(bucketIndex,seq,pos,firSeqFlag)
                firSeqFlag = False
                #print(dna.reverse_complement(seq))
                #print("finish a encode seq")
                self.addtoGraph(bucketIndex,seq,pos)
                #print("finish add a seq")
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
    def Successors_in_graph(self,bucketIndex,K):
    	succ = []
    	#return succ
    	if K in self.bucketDict:
	    	for base in range(4):
	    		if self.bucketDict[K][base] > 0 and self.sufficient(bucketIndex,K,base):
	    			succ.append(K[1:]+self.num2dna[base])
    	return succ
    def sufficient(self,bucketIndex,K,base):
        totalSum = sum(self.bucketDict[K])
        if self.bucketDict[K][base]/totalSum < self.deleteBifurRatio:
            return False
        return True
    def encodeSeq(self,bucketIndex,seq,pos, firSeqFlag):
        if firSeqFlag : # encode the first sequence in bucket
            for base in seq[:pos]:
                self.firstSeq.append('0b00' if base.upper() == 'A'  else  '0b01' if base.upper() == 'C' else '0b10' if base.upper() == 'G' else '0b11')
            for base in seq[pos+self.kmerLen:]:
                self.firstSeq.append('0b00' if base.upper() == 'A'  else  '0b01' if base.upper() == 'C' else '0b10' if base.upper() == 'G' else '0b11')
        else:
            #encode the left part
            encodePos = 0
            currentPos = 0
            K = dna.reverse_complement(seq[pos:self.kmerLen+pos])
            for base in dna.reverse_complement(seq[:pos]):
                K_next = self.Suffix(K) + base
                pred = self.Successors_in_graph(bucketIndex,K)
                if len(pred) == 1:
                    if pred[0] != K_next:
                        self.newNodeNum += 1
                        self.numFlag.append('0b1')
                        currentPos = encodePos
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
                        self.cutBifur(bucketIndex,K)
                        self.setFreqs(K,pred)
                        self.encodeSeqPathL.write(self.freqs,self.dna2num[base])
                    K = K_next
                encodePos += 1
            #encode the right part
            encodePos = 0
            currentPos = 0
            K = seq[(self.indexLen+pos-self.kmerLen):(self.indexLen+pos)]
            for base in seq[pos+self.indexLen:]:
                K_next = self.Suffix(K) + base
                succ = self.Successors_in_graph(bucketIndex,K)
                if len(succ) == 1:
                    if succ[0] != K_next:
                        self.newNodeNum += 1
                        self.numFlag.append('0b1')
                        currentPos = encodePos
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
                        self.cutBifur(bucketIndex,K)
                        self.setFreqs(K,succ)
                        self.encodeSeqPathR.write(self.freqs, self.dna2num[base])
                    K = K_next
                encodePos += 1 

        return
    def setFreqs(self,K,next):
        for i in range(4):
            self.freqs.set(i,1)
        for base in next:
            self.freqs.set(self.dna2num[base[-1]],int(self.bucketDict[K][self.dna2num[base[-1]]])+1)
    def cutBifur(self,bucketIndex,K):
        return
        totalSum = sum(self.bucketDict[K]) 
        for base in range(4):
            if (self.bucketDict[K][base]/totalSum) < self.deleteBifurRatio:
                self.bucketDict[K][base] = 0
        return
    def outputNodeInfo(self):
    	print("the numbers of simpleNodeNum is %d, take up %.2f%%"%(self.simpleNodeNum,self.simpleNodeNum*100.0/(self.simpleNodeNum+\
    		  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
    	print("the numbers of newNodeNum is %d, take up %.2f%%"%(self.newNodeNum,self.newNodeNum*100.0/(self.simpleNodeNum+\
    		  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
    	print("the numbers of tipNodeNum is %d, take up %.2f%%"%(self.tipNodeNum,self.tipNodeNum*100.0/(self.simpleNodeNum+\
    		  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
    	print("the numbers of bifurNodeNum is %d, take up %.2f%%"%(self.bifurNodeNum,self.bifurNodeNum*100.0/(self.simpleNodeNum+\
    		  self.newNodeNum + self.tipNodeNum + self.bifurNodeNum)))
    def addtoGraph(self,bucketIndex,seq,pos):
        for kmer in self.getKmerR(seq,pos):
            if str(kmer[:-1]) not in self.bucketDict:
              self.bucketDict.setdefault(str(kmer[:-1]), [0,0,0,0])
            self.bucketDict[str(kmer[:-1])][self.dna2num[str(kmer[-1])]] += 1
        for kmer in self.getKmerL(seq,pos):
            if dna.reverse_complement(kmer[1:]) not in self.bucketDict:
                self.bucketDict.setdefault(dna.reverse_complement(kmer[1:]), [0,0,0,0])           
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
        self.bitoutL.close()
        self.bitoutR.close()
        leftsize=os.path.getsize(self.outPutPath + ".bifurL")
        rightsize=os.path.getsize(self.outPutPath + ".bifurR")
        firstSeqsize = os.path.getsize(self.outPutPath + ".firSeq.lz")
        numFlagsize = os.path.getsize(self.outPutPath + ".numFlag.lz")
        filesize=[leftsize,rightsize,firstSeqsize,numFlagsize]
        print("[bifurL,bifurR,firSeq,numFlag]")
        print(filesize,"sum",sum(filesize),"bytes")
        return





