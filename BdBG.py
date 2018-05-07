from __future__ import print_function

import sys
import argparse
import os

from readFQ import get_chunks,get_chunksData
import copy
from bucket import encodeBucketClass, decodeBucketClass
from deBruijnGraph import encodeGraphClass, decodeGraphClass
def args_check(args):
    if not args.encode and not args.decode:
        sys.exit("you must give a -e or -d for encode/decode")
    if not args.input:
        sys.exit("you must give a file input with -i input")
    if not args.output:
        sys.exit("you must give a file output with -o output") 

def main(args):
    args_check(args)
    if args.encode:       
        sortseq = encodeBucketClass()
        sortseq.setBucketIndex(args.kmer)    
        sortseq.setKmerLen(args.kmer)
        sortseq.setBlockSize(args.block)
        sortseq.getReadLen(args.input)
        sortseq.setPath(args.input, args.output)
        sortseq.initialFile()
        graph = encodeGraphClass(args.kmer,args.kmer,args.output)
        graph.setOutPath(args.output)
        graph.initialOutFile()

        #sortseq.openInput(filename)
        for start, block in get_chunks(args.input,args.block):
            data = get_chunksData(args.input, start, block).split('\n')
            sortseq.setBuffer(data)
            for record in sortseq.fastq_iter():
                sortseq.replaceN()
                sortseq.sortSeqence()
                sortseq.recordNum += 1
                if args.verbose:
                    if sortseq.recordNum % 10000 == 0:
                        print('process %d records'% sortseq.recordNum)
            ######rebuild table for singleton read
            sortseq.reassigned()    
            sortseq.outPutSeqence()
            print("finish encodeBucket")
            if args.verbose:
                sortseq.outputInfo()
            sortseq.saveSeqTable()
            del sortseq.mutiDict
            graph.setSequenceTable(sortseq.sequenceTableSave)
            print('process %d records'% sortseq.recordNum)
            sortseq.emptyPara()
            graph.encodeBucket()
            graph.outputEncodePath()
            print("finish outputEncodePath")
            if args.verbose:
                graph.outputNodeInfo()
            graph.emptyPara()

    ###compress bucket meta data
        sortseq.compressFile()
        graph.compressFile()
        
        del sortseq
        del graph
    #decode
    else:
        decodeBucket = decodeBucketClass()
        decodeBucket.setBucketIndex(args.kmer)
        decodeBucket.setPath(args.input, args.output)
        decodeBucket.decompress()

        decodeGraph = decodeGraphClass(args.kmer,args.kmer,args.input, args.output)
        decodeGraph.setMetaInfor(decodeBucket.readLen, decodeBucket.bucketIndexLen)
        decodeGraph.loadBucktData(decodeBucket.bucketIndex, decodeBucket.bucketCov,\
                                  decodeBucket.readIndexPos,decodeBucket.readrc, decodeBucket.readN)
        decodeGraph.decompress()
        del decodeBucket
        del decodeGraph


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'BdBG')
    parser.add_argument("-e", "--encode",
                    help="encoding",action="store_true")
    parser.add_argument("-d", "--decode",
                    help="decoding",action="store_true")
    parser.add_argument("-i", "--input",type=str,
                    help="inputFile")
    parser.add_argument("-o", "--output",
                    help="outputFile")
    parser.add_argument("-k", "--kmer",type=int, default=15,
                    help="kmer size for bucket and de Bruijn graph, default=15")
    parser.add_argument("-b", "--block",type=int,default=1024,
                    help="input block size (M), default=1024")
    parser.add_argument("-v","--verbose", action="store_true",
                    help="verbose information")
    args = parser.parse_args()
    main(args)







