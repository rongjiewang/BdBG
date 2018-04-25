# BdBG
BdBG - BdBG: a bucket-based method for bacterial genome sequencing data compression with a dynamic de Bruijn graph.

Currently, BdBG works with FASTQ files and supports only compression of DNA stream, discarding the read names and qualities. By compressing only the DNA stream, it can reduce the bacterial genome sequencing data from 8 bits per base to 0.23 bits per base, achieving compression ratio 2.875%.


# Usage

DNA stream compression using BdBG is a 2 stage process, consisting of bucket and de Bruijn raph subprograms in chain. However, to decompress the DNA stream, only de Bruijn graph is needed.

## BdBG.py

BdBG.py first performs DNA records clustering into separate bucket representing signatures. As an input it takes a single or a set of FASTQ files and stores the output to five separate files: `*.index.lz` file, containing encoded bucktes index stream, `*.cov.lz` file, containing the number of reads in the buckets, `*.indexPos.lz` file, containing the bucket index positons in each read, `*.rc` file, containing wether the read in forward or backward, and `*.N.lz` file, containing the character N postions and length in the reads.

And then, BdBG.py performs encoding the read as a path in the dynamic de Bruijn graph in each bucket independently. It stores the output to four separate files: `*.bifurL` file, containing encode the read left bifurcation path from the beginning 'anchor' k-mer, `*.bifurR` file, containing encode the read right bifurcation path from the beginning 'anchor' k-mer, `*.firSeq.lz` file, containing the first reads in the buckets, and `*.numFlag.lz` file, containing the New nodes positions flag in the bifuraction list.

### Command line
BdBG.py is run from the command prompt:

    python BdBG.py <-e|-d> [options]

in one of the two modes:
* `-e` - encoding,
* `-d` - decoding,

with available options:
* `-i<file>` - input file,
* `-o<f>` - output files prefix,
* `-k<n>` - k-mer length, default: `15`,
* `-b<n>` - FASTQ input buffer size (in MB), default: `1024`,
* `-v` - verbose mode, default: `false`.


### Examples
Encode reads from `test.fastq` file, using signtature length of `15` and `1024` MB FASTQ block buffer, saving output to `output ` files:

    python BdBG.py -e -i test.fastq -o output 
    
Decode reads from `output` files and save the DNA reads to `input.dna` file.

    python BdBG.py -d -i output -o input.dna

