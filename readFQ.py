from __future__ import print_function

import sys
import csv
import screed
import argparse
import textwrap
import os
import string
import time
import numpy as np
from ngs_plumbing import dna ## DNA to 2bit
from bitstring import * ##glombe integer
from itertools import islice
import concurrent.futures

from khmer import __version__
from khmer.khmer_args import (sanitize_help, ComboFormatter, info,
                              _VersionStdErrAction)
import threading
from Bio import SeqIO

mutex = threading.RLock()

def readfq(filename): # this is a generator function
    fp = open(filename,'r')
    start = fp.tell()
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    #if (fp.tell()-start) >= blockSize:
                        #print("out side, tell %d, start %d,blockSize %d "%(fp.tell(),start,blockSize))
                        #break
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

class Reader(threading.Thread):
    def __init__(self, num):
        super(Reader,self).__init__()
        self.num = num

    # def run(self):
    #     size = mb * (1024**2)
    #     f = open(file,'r')
    #     while True:
    #         with mutex:
    #             next_n_lines = list(islice(fp, 4))
    #             if not next_n_lines:
    #                 return
    #             print('%d:%s' % (self.num, line), end='')

def get_chunks(file, mb=256):
	size = mb * (1024**2)
	f = open(file,'r')
	while 1:
		with mutex:
			start = f.tell()
			f.seek(size, 1)
			line = f.readline()
			if not line:
				yield start,f.tell()-start
				break 
			while line and not line.startswith('@'):
				line = f.readline()
			line2 = f.readline()
			if not line2:
				yield start,f.tell()-start
				break 
			if not line2.startswith('@'):
				f.seek(-len(line2),1)
			else:
				line = line2
			#else:
			# back up the length of the fasta header
			f.seek(-len(line), 1)
			# tuple up
			yield start,f.tell()-start
	f.close()

def parse(y):
	for x in range(0,(len(y)-1),4):
		yield y[x],y[x+1],y[x+2]
def getRecords(y):
	records = []
	for x in range(0,(len(y)-1),4):
		records.append([y[x],y[x+1],y[x+2]])
	return records
def main(filename):
	recordNum = 0
	datasize = 0
	for start, block in get_chunks(filename):
		#data = readchunks(filename,a,b)
		#out.write(data)
		#print(start, block)
		data = get_chunksData(filename, start, block)
		#print ("data is %s"%data)
		for name,seq,qual in parse(data):
			# print(name)
			# print(seq)
			# print(qual)
			recordNum += 1
			if(recordNum % 10000 == 0):
				print(recordNum)
		datasize += len(data)
	print(threading.currentThread().getName())
	print('data size is %d, record num is %d'%(datasize, recordNum))

def get_chunksData(file,offset,blockSize):
	f = open(file,'r')
	f.seek(offset)
	return f.read(blockSize)
	f.close()
def process(data):
	return data

if __name__ == '__main__':
	filename = "/home/wang/rjwang/data/benchmark/FQ/1SRR554369_1.fastq"
	#filename = "/home/rjwang/e/soft/mypthon/mince/input/mda-test.fastq"
	#filename = "/home/rjwang/compress/Data/read/Pseudomonas/SRR554369_1.fastq"
	#filename = "/home/rjwang/e/data/benchmark/FQ/5ERR174310_1.fastq"



	# datasize = 0
	# recordNum = 0
	# for start, block in get_chunks(filename):
	# 	#data = readchunks(filename,a,b)
	# 	#out.write(data)
	# 	#print(start, block)
	# 	data = get_chunksData(filename, start, block)
	# 	#print ("data is %s"%data)
	# 	for name,seq,qual in parse(data):
	# 		# print(name)
	# 		# print(seq)
	# 		# print(qual)
	# 		recordNum += 1
	# 		if(recordNum % 100000 == 0):
	# 			print(recordNum)
	# 	datasize += len(data)
	# print('data size is %d, record num is %d'%(datasize, recordNum))


	# recordnum = 0
	# for data in readfq(filename):
	# 	recordnum +=1
	# 		if(recordNum % 100000 == 0):
	# 	print(recordNum)
	# print("records num is %d"%(recordnum))#24S






	# input_iter = screed.open(filename)
	# recordNum = 0
	# for record in input_iter:
	# 	recordNum += 1
	# print("process %d records"% recordNum)###152S


	# time.sleep(3)
	# recordNum = 0
	# with open(filename) as f:
	# 	while True:
	# 		next_n_lines = list(islice(f, 4000000))
	# 		if not next_n_lines:
	# 		    break
	# 		# process next_n_lines
	# 		for line in next_n_lines:
	# 			recordNum += 1
	# 			print(line)
	# 		if recordNum % 4000000 == 0:
	# 			print(time.ctime())
	# 			x = recordNum / 4
	# 			print("process %d records"% x)
	# fp = open(filename,'r')
	# thread = 1
	# t=[]
	# for i in range(thread):
	# 	my_thread = threading.Thread(target=main, args=(fp,filename,))
	# 	t.append(my_thread)
	# 	t[i].start()
	# for i in range(thread):
	# 	t[i].join()



	datasize = 0
	recordNum = 0
	for start, block in get_chunks(filename):
		#data = readchunks(filename,a,b)
		#out.write(data)
		#print(start, block)
		data = get_chunksData(filename, start, block).split('\n')
		#records = getRecords(data)
		datasize += len(data)
		#del data
		#print ("data is %s"%data)
		for name,seq,qual in parse(data):
			# print(name)
			# print(seq)
			# print(qual)
			#process(seq)
			recordNum += 1
			if(recordNum % 100000 == 0):
				print(recordNum)
		# with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
		# 	for seq, result in zip(data, executor.map(process, data)):
		# 		#print('%s is prime: %s' % (seq, result))
		# 		recordNum += 1
		# 		if(recordNum % 100000 == 0):
		# 			print(recordNum)
		
	print('data size is %d, record num is %d'%(datasize, recordNum))#10.5S

# def read_in_chunks(file_object, chunk_size=1024):
#     while True:
#         data = file_object.read(chunk_size)
#         if not data:
#             break
#         yield data


# 	f = open(file, 'rb')
# 	for piece in read_in_chunks(f):
# 	    process_data(piece)          
# 	f.close()

	# recorderNum = 0
	# for seq_record in SeqIO.parse(filename, "fastq"):
	# 	recorderNum +=1
 # 		if(recordNum % 100000 == 0):
 # 			print(recordNum)
	# print (recorderNum)### 192S