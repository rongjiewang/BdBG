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

import sys
import argparse

from bucket import encodeBucketClass, decodeBucketClass
from deBruijnGraph import encodeGraphClass, decodeGraphClass

def args_check(args):
    if not args.encode and not args.decode:
        sys.exit("you must give a -e or -d for encode/decode")
    if not args.input and not args.paired:
        sys.exit("you must give a file input with -i input for single end data or -p -1 input1 -2 input2 for paired-end data")
    if not args.output:
        sys.exit("you must give a file output with -o output")
    return 

def main(args):
    args_check(args)
    #encode
    if args.encode:       
        en_bucket = encodeBucketClass(args.input, args.output, args.paired, \
                                        args.input1, args.input2, args.kmer, args.lossless, args.verbose)
        en_bucket.encode()

        en_graph = encodeGraphClass(args.output, args.paired, args.kmer, \
                                        args.verbose, en_bucket.sequenceTableSave)
        del en_bucket

        en_graph.encode()

        del en_graph

        sys.exit()

    #decode
    else:
        de_bucket = decodeBucketClass(args.input, args.output, args.verbose)

        de_bucket.decode()

        de_graph = decodeGraphClass(args.input, args.output, de_bucket.paired, de_bucket.readNum,\
                                    de_bucket.bucketIndexLen, de_bucket.lossless, de_bucket.verbose)

        de_graph.loadBucktData(de_bucket.bucketIndex, de_bucket.bucketCov, de_bucket.readIndexPos,\
                                de_bucket.readrc, de_bucket.readN, de_bucket.readLen, de_bucket.readOrder)
        del de_bucket

        de_graph.decode()

        del de_graph

        sys.exit()

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
    parser.add_argument("-p", "--paired",
                    help="paired-end flag",action="store_true")
    parser.add_argument("-1", "--input1",
                    help="paired-end file1")
    parser.add_argument("-2", "--input2",
                    help="paired-end file2")
    parser.add_argument("-l", "--lossless",
                    help="keep the reads orders, default:false, \
                    if encode paired-end files, default:ture ",action="store_true")
    parser.add_argument("-k", "--kmer",type=int, default=15,
                    help="kmer size for bucket and de Bruijn graph, default=15")
    parser.add_argument("-v","--verbose", action="store_true",
                    help="verbose information")

    args = parser.parse_args()

    main(args)







