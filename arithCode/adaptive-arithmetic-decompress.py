# 
# Decompression application using adaptive arithmetic coding
# 
# Usage: python adaptive-arithmetic-decompress.py InputFile OutputFile
# This decompresses files generated by the adaptive-arithmetic-compress.py application.
# 
# Copyright (c) Project Nayuki
# 
# https://www.nayuki.io/page/reference-arithmetic-coding
# https://github.com/nayuki/Reference-arithmetic-coding
# 

import sys
import arithmeticcoding
python3 = sys.version_info.major >= 3


# Command line main application function.
def main(args):
	# Handle command line arguments
	if len(args) != 2:
		sys.exit("Usage: python adaptive-arithmetic-decompress.py InputFile OutputFile")
	inputfile  = args[0]
	outputfile = args[1]
	
	# Perform file decompression
	bitin = arithmeticcoding.BitInputStream(open(inputfile, "rb"))
	with open(outputfile, "wb") as out:
		try:
			decompress(bitin, out)
		finally:
			bitin.close()


def decompress(bitin, out):
	initfreqs = arithmeticcoding.FlatFrequencyTable(257)
	freqs = arithmeticcoding.SimpleFrequencyTable(initfreqs)
	dec = arithmeticcoding.ArithmeticDecoder(bitin)
	while True:
		# Decode and write one byte
		symbol = dec.read(freqs)
		if symbol == 256:  # EOF symbol
			break
		out.write(bytes((symbol,)) if python3 else chr(symbol))
		freqs.increment(symbol)


# Main launcher
if __name__ == "__main__":
	main(sys.argv[1 : ])
