# BdBG
BdBG - BdBG: a bucket-based method for bacterial genome sequencing data compression with a dynamic de Bruijn graph.

Currently, BdBG works with FASTQ files and supports only compression of DNA stream, discarding the read names and qualities. By compressing only the DNA stream, it can reduce the bacterial genome sequencing data from 8 bits per base to 0.23 bits per base, achieving compression ratio 2.875%.

## Install
This is a step by step instruction for installing the BdBG on linux.
### Requirements for python 2.7* modules
* Screed
* numpy
* bitstring

### Command to install
    pip install -r requirements.txt

## Usage

DNA stream compression using BdBG is a 2 stage process, consisting of bucket and de Bruijn raph subprograms in chain. However, to decompress the DNA stream, only de Bruijn graph is needed.

### BdBG.py

BdBG.py first performs DNA records clustering into separate bucket representing signatures. As an input it takes a single or a set of FASTQ files and stores the output to five separate files: `*.index.lz` file, containing encoded bucktes index stream, `*.cov.lz` file, containing the number of reads in the buckets, `*.indexPos.lz` file, containing the bucket index positons in each read, `*.rc` file, containing wether the read in forward or backward, and `*.N.lz` file, containing the character N postions and length in the reads.

And then, BdBG.py performs encoding the read as a path in the dynamic de Bruijn graph in each bucket independently. It stores the output to four separate files: `*.bifurL` file, containing encode the read left bifurcation path from the beginning 'anchor' k-mer, `*.bifurR` file, containing encode the read right bifurcation path from the beginning 'anchor' k-mer, `*.firSeq.lz` file, containing the first reads in the buckets, and `*.numFlag.lz` file, containing the New nodes positions flag in the bifuraction list.

There are a option parameter '-l' to unchange the raw reads orders by extra adding a `*.order.lz` file to compressed files. The parameter '-l' default is false. However, when encode paired-end files, default value become ture, beceause we have to keep the paired information for paired-end reads. 


### Command line
BdBG.py is run from the command prompt:

    python BdBG.py <-e|-d> [options]

in one of the two modes:
* `-e` - encoding,
* `-d` - decoding,

with available options:
* `-i<file>` - input.fastq file,
* `-o<f>` - output files prefix,
* `-p` - paired-end file flag,
* `-1<file>` - input_1.fastq file,
* `-2<file>` - input_2.fastq file,
* `-l` - lossless encode the read, means keep the reads order, default:false, if encode paired-end files, default:ture,
* `-k<n>` - k-mer size, default: `15`,
* `-v` - verbose mode, default: `false`.


### Examples
Encode single-end reads with `test.fastq` file, output with prefix `encode_test`:

    python BdBG.py -e -i test.fastq -o encode_test
    
Encode single-end reads with `test.fastq` file and keep the raw read orders, output with prefix `encode_test`:

    python BdBG.py -e -l -i test.fastq -o encode_test
    
Encode paired-end reads with `test_1.fastq` file and `test_2.fastq` file, output with prefix `encode_test`:

    python BdBG.py -e -p -1 test_1.fastq -2 test_2.fastq -o encode_test 
    
Decode reads from output prefix `encode_test` files and save the DNA reads with prefix `decode_test`.

    python BdBG.py -d -i encode_test -o decod_test
    
### Verify the compress & decompress correctness
For verify the correctness, first extarct dna sequence from fastq file, by using the file ./scripts/extract_dna_from_fastq.sh.

    sh ./scripts/extract_dna_from_fastq.sh input.fastq input.dna 
    
Then, sort the input.dna and output.dna. Last, check the difference beween input.dna and output.dna.

    sort input.dna > sorted_input.dna
    sort ouput.dna > sorted_output.dna
    diff sorted_input.dna sorted_output.dna
    
If compress the raw sequence with parameter '-l', which means to unchange the reads order, you can compare the input.dna and output.dna directly.

    diff input.dna output.dna
    
If there nothing out in creeen with shell command diff, proof that the compresstion & decompression process is correct.
    
## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
The code of BdBG was based in part on the source code of the Arithmetic package [Nayuki, 2016](https://github.com/nayuki/Reference-arithmetic-coding).

## Contact
If you have any question, please contact the author rjwang.hit@gmail.com
