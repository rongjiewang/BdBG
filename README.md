# BdBG
BdBG - BdBG: a bucket-based method for compressing genome sequencing data with dynamic de Bruijn graphs.

BdBG, a new alignment-free and reference-free compression of FASTQ sequences stream, the method is based on the concept of bucketing similar reads into the same bucket and compressing reads in each bucket individually by a dynamic de Bruijn graph. Experimental results on eight different genome and transcriptome sequencing datasets testified the compression performance of BdBG is better than that of GZIP, LEON and MINCE, with improvements of up to 83%, 81%, and 52%, respectively.

## Downloading datasets
### Dataset ID  &   direct link
    ERR1147042	ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/ERX/ERX122/ERX1225844/ERR1147042
    ERR034088	ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/ERX/ERX012/ERX012615/ERR034088
    SRR554369	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR554/SRR554369
    SRR959239	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR959/SRR959239
    ERR418881	ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/ERX/ERX385/ERX385178/ERR418881
    MH0001.081026	http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.1.fq.gz
                    http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.2.fq.gz
    SRR327342	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR327/SRR327342
    SRR037452	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR037/SRR037452
    
## Install
This is a step by step instruction for installing the BdBG for python 2.7*.
### Requirements for python modules & versions
* screed >= 1.0
* bitstring >= 3.1.5
* bitstream >= 2.5.4
* bitio >= 0.2
### Requirements for text compression tool & version
* plzip >= 0.9

### Command to install
    pip install -r requirements.txt
    sudo apt-get install plzip

## Usage

Sequence data stream compression using BdBG is a 2 stage process, consisting of bucket and de Bruijn raph subprograms in chain. However, to decompress the DNA stream, only de Bruijn graph is needed.

### BdBG OutPut

BdBG first performs reading bucket, output streams contains five files: 

* `*.index.lz`： bucktes index stream；
* `*.cov.lz`： the number of reads in the buckets;
* `*.indexPos.lz`: the bucket index positons in each read;
* `*.rc`:  whether the read in forward or in reverse-complement direction;
* `*.N.lz`: the characters "N" postions and length in the reads.

Then, BdBG encodes the read as a path in the dynamic de Bruijn graph in each bucket independently. output streams contains five files: 

* `*.bifurL`: the read left bifurcation path from the beginning 'anchor' k-mer;
* `*.bifurR`: the read right bifurcation path from the beginning 'anchor' k-mer;
* `*.firSeq.lz`: the first read in each bucket; 
* `*.numFlag.lz`: the New node position flags in the bifuraction list.
* `*.order.lz`: reserve the raw read orders for paired-end reads. it is a option parameter '-l' for single-end reads, defaut:false. 


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
