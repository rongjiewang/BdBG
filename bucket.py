"""MIT License

Copyright (c) 2018 rongjiewang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

from __future__ import print_function
import sys
import screed
import os
import string
from bitstring import * #for BitStream()
from bitstream import BitStream as sream

from screed import Record
from utils import to_str
import copy
import struct
import math
import dna #for reverse_complement

class encodeBucketClass(object):
    """encoding for bucket sequence"""
    def __init__(self, input, output, ispaired, input1, input2, kmerLen, lossless, verbose):
        self.mutiDict = {} ### kmers for buckets 
        self.bucketTable = {} ### buckets index and it's reads number
        self.sequenceTable = {} ####bucket as diction index, within bucket include (seq,reverseFlag,indexPos) 
        self.sequenceTableSave = {}
        self.read = {'sequence':"", 'reverse':"", 'indexPos':"", 'N':"", 'len':"", 'order':""}
        self.kmerLen = kmerLen
        self.bucketIndexLen = kmerLen
        # self.blockSize = 1024
        self.recordNum = 0
        self.onsies = 0 ####the number of bucket which only have one seq
        self.leftOnesies = 0 ## the number onesies after re-assigned
        self.sequenceN = []  ### the info of N in read  present by (readID,pos,len)
        self.reversedNum = 0 ##total reverse number
        self.skipZone = 0
        self.record = Record() #read record
        self.buffer = None
        self.sequenceNFlag = False #indicate if the seq include N
        self.input = input
        self.fileOutputPath = output
        self.seqLen = 0
        self.fileoutName = {"index":"", "cov":"", "rc":"", "N":"", "indexPos":"", "len":0, "order":0}
        self.paired = ispaired
        self.input1 = input1
        self.input2 = input2
        self.verbose = verbose
        self.maxKmer = "T"*kmerLen
        self.lossless = lossless
        self.removeOutputFile()

    def setskipZone(self,skip):
        self.skipZone = skip
        return

    def getReadLen(self,filename):
        input_iter = screed.open(filename)
        for record in input_iter:
            self.seqLen = len(record.sequence)
            break
        return

    def mergePairRead(self, record1, record2):
        self.record['sequence'] = str(record1['sequence']) + dna.reverse_complement(str(record2['sequence']))
        return

    def setBucketIndex(self,bucketIndexLen):
        self.bucketIndexLen = bucketIndexLen
        return

    def getBucketKmers(self,seq):
        for i in xrange(self.seqLen-self.skipZone-self.bucketIndexLen+1):
            yield seq[i:i+self.bucketIndexLen]

    def getReadKmers(self,seq):
        for i in xrange(self.seqLen-self.kmerLen+1):
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
        try:
            self.mutiDict[index]
            for kmer in self.getReadKmers(read):
                try:
                    self.mutiDict[index][kmer]
                    overlapNum += 1
                except KeyError:
                    continue
            return overlapNum
        except KeyError:
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

    def addTomutiDict(self,index,read):
        try:
            self.mutiDict[index]
            for kmer in self.getReadKmers(read):
                self.mutiDict[index][kmer] = self.mutiDict[index].get(kmer,0) + 1
        except KeyError:
            self.mutiDict.setdefault(index, {})
            for kmer in self.getReadKmers(read):
                self.mutiDict[index][kmer] = self.mutiDict[index].get(kmer,0) + 1
        return

    def getMostCommenIndex(self):
        read = self.record['sequence']
        rcread = dna.reverse_complement(str(read))
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
        rcread = dna.reverse_complement(str(read))
        minKmer = self.maxKmer
        reverseFlag = False
        indexPos = 0
        kmerIter = 0
        for kmer in self.getBucketKmers(read):
            if kmer < minKmer:
                minKmer = kmer
                indexPos = kmerIter
            kmerIter +=1

        kmerIter = 0
        for kmer in self.getBucketKmers(rcread):
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
            self.read['sequence'] = dna.reverse_complement(str(self.record['sequence']))
        else:
            self.read['sequence'] = str(self.record['sequence'])
        self.read['reverse'] = reverseFlag
        self.read['indexPos'] = indexPos
        self.read['len'] = len(self.record['sequence'])
        self.read['order'] = self.recordNum
        self.sequenceTable.setdefault(index, []).append(copy.deepcopy(self.read))
        return

    def reSortSeqence(self):##for reasin single read 
        """Run over the given file and sort the sequences."""
        index,reverseFlag,indexPos = self.getMostCommenIndex()
        self.bucketTable[index] = self.bucketTable.get(index,0) + 1
        if reverseFlag:
            self.read['sequence'] = dna.reverse_complement(str(self.record['sequence']))
        else:
            self.read['sequence'] = str(self.record['sequence'])
        self.read['reverse'] = reverseFlag
        self.read['indexPos'] = indexPos
        self.sequenceTable.setdefault(index, []).append(copy.deepcopy(self.read))
        return

    def reassigned(self): 
        for countIndex in self.bucketTable.keys():
            if self.bucketTable[countIndex] == 1:
                if self.sequenceTable[countIndex][0]['reverse']:##recover the raw sequence
                    self.record['sequence'] = dna.reverse_complement(str(self.sequenceTable[countIndex][0]['sequence']))
                else:
                    self.record['sequence'] = self.sequenceTable[countIndex][0]['sequence']
                self.read['N'] =  self.sequenceTable[countIndex][0]['N']## take N infor for reassign
                self.read['order'] = self.sequenceTable[countIndex][0]['order']
                self.read['len'] = self.sequenceTable[countIndex][0]['len']
                self.seqLen = len(self.record['sequence'])
                del self.mutiDict[countIndex]
                del self.bucketTable[countIndex]
                del self.sequenceTable[countIndex]
                self.reSortSeqence()
        return

    def outputInfo(self):
        for countIndex in self.bucketTable.keys():
            if self.bucketTable[countIndex] == 1:
                self.onsies += 1
        print("the number onsies is",self.onsies)
        print("the number of bucket is",len(self.bucketTable.keys()))
        print("the number of read is",self.recordNum)
        print("the number of reversed reads is {0},occupy {1}% ".format(self.reversedNum, self.reversedNum*100/self.recordNum))
        return

    def sort_key(self,seq):
        return seq['indexPos'],seq['sequence'][(seq['indexPos']+self.kmerLen):]

    def setKmerLen(self,L):
        self.kmerLen = L
        return

    def removeOutputFile(self):
        for name in self.fileoutName:
            file = self.fileOutputPath + "." + name                            
            if os.path.exists(file):
                os.remove(file)
        return

    def outPutSeqence(self):
        for name in self.fileoutName.keys()[:-1]:
            file = self.fileOutputPath + "." + name                              
            self.fileoutName[name] = open(file,"a+")
        if self.lossless:
            name = self.fileoutName.keys()[-1]
            file = self.fileOutputPath + "." + name
            self.fileoutName[name] = open(file,"a+")

        ##store the meta information about read isPaired (1 byte), bucket index length (1 byte) and reads number (4 byte).
        self.fileoutName["index"].write(chr(self.paired))
        self.fileoutName["index"].write(chr(self.bucketIndexLen))
        self.fileoutName["index"].write(chr(self.recordNum/16777216))
        self.fileoutName["index"].write(chr(self.recordNum%16777216/65536))
        self.fileoutName["index"].write(chr(self.recordNum%16777216%65536/256))
        self.fileoutName["index"].write(chr(self.recordNum%256))
        self.fileoutName["index"].write(chr(self.lossless))#indicate whether encode the read in keep the read order way

        ###output bucket########
        curren_index = 0
        for index in sorted(self.bucketTable):
            indexBin = self.getIndexMinKmer(index) -curren_index
            curren_index = indexBin + curren_index
            self.fileoutName["index"].write(str(indexBin)+"\n")###store index in number
            self.fileoutName["cov"].write(str(self.bucketTable[index]-1)+"\n")###store covrage in number
        self.fileoutName["index"].close()
        self.fileoutName["cov"].close()

        ###output sequence information#######
        RCbin = BitStream()
        Nindex = BitStream()
        for index in sorted(self.sequenceTable.keys()):
            prePos = 0  ###for delta encode index position
            for read in sorted(self.sequenceTable[index],key=self.sort_key):
                if read['reverse']:
                    self.reversedNum += 1
                    RCbin.append('0b1')
                else:
                    RCbin.append('0b0')
                self.fileoutName["indexPos"].write(chr(read['indexPos']-prePos))
                prePos = read['indexPos']
                if read['N']:
                    Nindex.append('0b1')
                else:
                    Nindex.append('0b0')
                self.fileoutName["len"].write(chr(read['len'])) # save read length
                if self.lossless:
                    self.fileoutName["order"].write(struct.pack('I', read['order']))## save read order infor

        RCbin.tofile(self.fileoutName["rc"]) 
        self.fileoutName["rc"].close()
        self.fileoutName["indexPos"].close()
        self.fileoutName["len"].close()
        if self.lossless:
            self.fileoutName["order"].close()
        
        ####output N pos and Length information#######
        Nindex.tofile(self.fileoutName["N"])
        for index in sorted(self.sequenceTable.keys()):
            for read in sorted(self.sequenceTable[index],key=self.sort_key):
                if read['N']:
                    self.fileoutName["N"].write(chr(len(read['N'][1]))) ## N numbers
                    prePos = 0
                    for pos in read['N'][1]:
                        self.fileoutName["N"].write(chr(pos-prePos))
                        prePos = pos
                    for posLen in read['N'][2]:
                        self.fileoutName["N"].write(chr(posLen))
        self.fileoutName["N"].close()
        return

    def getNInRead(self):
        read = self.record['sequence']
        pos = []
        posLen = []
        i = 0        
        while(i < self.seqLen):
            if read[i] == "N":
                nLen = 1
                pos.append(i)
                i += 1
                while (i < self.seqLen ):
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
            self.read['N'] = False
        else:
            pos,length = self.getNInRead()
            self.read['N'] = (True, pos, length)
            self.record['sequence'] = string.replace(self.record['sequence'], "N", "C")
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
            iterseq = iter(sorted(self.sequenceTable[bucketIndex],key=self.sort_key))     
            for read in sorted(self.sequenceTable[bucketIndex],key=self.sort_key):
                self.sequenceTableSave.setdefault(bucketIndex, []).append([str(read['sequence']),read['indexPos']])
        return

    def compressFile(self):
        if self.lossless:
            for name in self.fileoutName:
                file = self.fileOutputPath + "." + name    
                if os.path.exists(file + ".lz"):
                    os.remove(file + ".lz")
                comand= "plzip " + file
                os.system(comand)
        else:
            for name in self.fileoutName.keys()[:-1]:
                file = self.fileOutputPath + "." + name    
                if os.path.exists(file + ".lz"):
                    os.remove(file + ".lz")
                comand= "plzip " + file
                os.system(comand)
        return

    def analysisFileSize(self):
        for name in self.fileoutName:
            file = self.fileOutputPath + "." + name    
            if os.path.exists(file + ".lz"):
                filesize = os.path.getsize(file + ".lz")
                print("file %s size is %d", (name, filesize))
        return

    def encode(self):
        # encode paired-end data
        if self.paired:
            self.lossless = True
            input_iter1 = screed.open(self.input1)
            input_iter2 = screed.open(self.input2)
            for record1, record2 in zip(input_iter1, input_iter2):
                self.record = record1
                self.seqLen = len(self.record['sequence'])
                if self.seqLen > 256:
                    sys.exit("BdBG compression tool for genomic reads. \
                        it works on fixed or variable length reads less than 256.")

                self.replaceN()
                self.sortSeqence()
                self.recordNum += 1
                if self.verbose:
                    if self.recordNum % 100000 == 0:
                        print('process %d records'% self.recordNum)
                self.record = record2
                self.seqLen = len(self.record['sequence'])
                if self.seqLen > 256:
                    sys.exit("BdBG compression tool for genomic reads. \
                        it works on fixed or variable length reads less than 256.")
                self.replaceN()
                self.sortSeqence()
                self.recordNum += 1
                if self.verbose:
                    if self.recordNum % 100000 == 0:
                        print('process %d records'% self.recordNum)
            ######rebuild table for singleton read
            self.reassigned()    
            self.outPutSeqence()
            if self.verbose:
                self.outputInfo()
            self.saveSeqTable()
            del self.mutiDict
            self.compressFile()
            print("finish bucket")            
            print('process %d records'% self.recordNum)

        # encode single-end data
        else:
            input_iter = screed.open(self.input)
            for record in input_iter:
                self.record = record
                self.seqLen = len(self.record['sequence'])
                if self.seqLen > 256:
                    sys.exit("BdBG compression tool for genomic reads. \
                        it works on fixed or variable length reads less than 256.")
                self.replaceN()
                self.sortSeqence()
                self.recordNum += 1
                if self.verbose:
                    if self.recordNum % 10000 == 0:
                        print('process %d records'% self.recordNum)
            ######rebuild table for singleton read
            self.reassigned()    
            self.outPutSeqence()
            if self.verbose:
                self.outputInfo()
            self.saveSeqTable()
            del self.mutiDict
            self.compressFile()
            print("finish bucket")            
            print('process %d records'% self.recordNum)
        return

class decodeBucketClass(object):
    """decoding for bucket sequence"""
    def __init__(self, input, output, verbose):
        self.bucketIndex = [] #bucket index 
        self.bucketCov = [] # reads number in bucket
        self.readIndexPos = [] #index positions in each read
        self.readOrder = [] # order positions for each raw read 
        self.readrc = sream() # read in forward or backward
        self.readN = {"flag":sream(), "pos":[], "l":[]} # N in read indicate, position and length
        self.bucketIndexLen = 0
        self.recordNum = 0
        self.fileInputPath = input
        self.fileOutputPath = output
        self.readLen = [] # set read length
        self.seqLen = 0 # current read length
        self.paired = False
        self.fileInputName = {"index":"", "cov":"", "rc":"", "N":"", "indexPos":"", "len":"", "order":""}
        self.num2dna = {0:"A",1:"C",2:"G",3:"T"}
        self.readNum = 0
        self.lossless = False
        self.verbose = verbose

    def setBucketIndex(self,bucketIndexLen):
        self.bucketIndexLen = bucketIndexLen
        return

    def setPath(self,inputPath,outputPath):
        self.fileInputPath = inputPath
        self.fileOutputPath = outputPath
        return

    def decode(self):
        self.unPlzipFiles()# unzip the meta file
        self.decodeIndex() # decode bucket index
        self.decodeCov() #decode the coverage for each bucket
        self.decodeIndexPos() #decode the index position for each read
        self.decodeRC() #decode the direction for each read
        self.decodeN() #decode the character "N" information for each read
        self.decodeLen() #decode the length for each read
        if self.lossless:
            self.decodeOrder() #decode the order information for each read
        self.deleteUnzipFiles() #delete unziped files
        return

    def unPlzipFiles(self):
        for name in self.fileInputName:
            file = self.fileInputPath + "." + name + ".lz"                           
            if os.path.exists(file):
                comand= "plzip -d -k " + file
                os.system(comand)
        return

    def decodeIndex(self):
        file = self.fileInputPath + ".index"
        indexNum = 0
        f = open(file,'r')
        # decode meta information about read length (1 byte), bucket index length (1 byte) and block size (M) (2 byte).
        self.paired = ord(f.read(1))
        self.bucketIndexLen = ord(f.read(1))
        self.readNum = ord(f.read(1))*16777216 + ord(f.read(1))*65536 + ord(f.read(1))*256 + ord(f.read(1))
        self.lossless = ord(f.read(1))
        for line in f.readlines():
            indexNum += int(line)
            self.bucketIndex.append(''.join(self.num2Base(indexNum)))
        f.close()
        return

    def num2Base(self, indexNum):
        indexBase = []
        while indexNum:
            indexBase.append(self.num2dna[(indexNum & 3)])
            indexNum = indexNum >> 2
        while len(indexBase) < self.bucketIndexLen:
            indexBase.append("A")
        return indexBase[::-1]

    def decodeCov(self):
        file = self.fileInputPath + ".cov"
        with open(file,'r') as f:
            for line in f.readlines():
                readNum = int(line)
                self.bucketCov.append(readNum)
        f.close()
        return

    def decodeIndexPos(self):
        file = self.fileInputPath + ".indexPos"
        with open(file,'r') as f:
            while 1:
                byte_s = f.read(1)
                if not byte_s:
                    break
                indexPos = ord(byte_s)
                self.readIndexPos.append(indexPos)
        f.close()
        return
    def decodeLen(self):
        file = self.fileInputPath + ".len"
        with open(file,'r') as f:
            while 1:
                byte_s = f.read(1)
                if not byte_s:
                    break
                length = int(ord(byte_s))
                self.readLen.append(length)
        f.close()
        return

    def decodeRC(self):
        file = self.fileInputPath + ".rc"
        with open(file,'rb') as f:
            data=f.read()
            self.readrc.write(data,bytes)
        f.close()
        return

    def decodeN(self):
        file = self.fileInputPath + ".N"
        with open(file,'rb') as f:
            Nflag=f.read(int(math.ceil(self.readNum/8.0)))
            self.readN["flag"].write(Nflag,bytes)## weather the read contain "N"
            byte_s = f.read(1) ### the number potions "N" the read have.
            while byte_s != "":
                pos = []
                length = []
                for i in range(ord(byte_s)):
                    pos.append(ord(f.read(1)))
                self.readN["pos"].append(pos)
                for i in range(ord(byte_s)):
                    length.append(ord(f.read(1)))
                self.readN["l"].append(length)
                byte_s = f.read(1)
        f.close()
        return

    def decodeOrder(self):
        file = self.fileInputPath + ".order"
        with open(file,'r') as f:
            while 1:
                byte_s = f.read(4)
                if not byte_s:
                    break
                order = struct.unpack('I', byte_s)
                self.readOrder.append(int(order[0]))
        f.close()
        return

    def deleteUnzipFiles(self):
        for name in self.fileInputName:
            file = self.fileInputPath + "." + name                           
            if os.path.exists(file):
                comand= "rm " + file
                os.system(comand)
        return

