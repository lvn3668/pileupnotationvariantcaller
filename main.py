# This is a sample Python script.
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# from typing import List

import argparse
from typing import Any, Dict

import numpy as np

import md5checksumgeneratorfordnabases.md5checksumgeneratorfordnasequences as md5
import pileupfilereader.pileupnotationreader as pileup
import referencegenomecheck.referencegenomecheck as refgenome
import chromosomelengthreader.readassemblyinfo as chrlengthreader
import samtoolspileupchecker.samtoolspileupchecker as smtools

# sequqence id , 1 based coordinate, reference base, number of reads covering that base,
# , / . matching positive / negative strand
#  ATGCN mismatch on fwd strand
# atgcn mismatch on rev strand
#

# If this is the first position covered by the read,
# a “^” character followed by the alignment's mapping quality encoded as an ASCII character.
# A single character indicating the read base and the strand to which the read has been mapped:
# Forward	Reverse	Meaning
# . dot	, comma	Base matches the reference base
# Base is a mismatch to the reference base
# >	<	Reference skip (due to CIGAR “N”) -> '>' on fwd strand and '<' on reverse strand
# *	*/#	Deletion of the reference base (CIGAR “D”) -> * deletion on fwd strand and # on reverse strand
# Deleted bases are shown as “*” on both strands unless --reverse-del is used,
# in which case they are shown as “#” on the reverse strand.
#
# If there is an insertion after this read base, text matching “\\+[0-9]+[ACGTNacgtn*#]+”:
# a “+” character followed by an integer giving the length of the insertion and then the inserted sequence.
# Pads are shown as “*” unless --reverse-del is used, in which case pads on the reverse strand will be shown as “#”.
# If there is a deletion after this read base, text matching “-[0-9]+[ACGTNacgtn]+”:
# a “-” character followed by the deleted reference bases represented similarly.
# (Subsequent pileup lines will contain “*” for this read indicating the deleted bases.)
# If this is the last position covered by the read, a “$” character.

# seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
# seq1 273 T 23  ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
# seq1 274 T 23  ,.$....,,.,.,...,,,.,...    7<7;<;<<<<<<<<<=<;<;<<6
# seq1 275 A 23  ,$....,,.,.,...,,,.,...^l.  <+;9*<<<<<<<<<=<<:;<<<<
# seq1 276 G 22  ...T,,.,.,...,,,.,....  33;+<<7=7<<7<&<<1;<<6<
# seq1 277 T 22  ....,,.,.,.C.,,,.,..G.  +7<;<<<<<<<&<=<<:;<<&<
# seq1 278 G 23  ....,,.,.,...,,,.,....^k.   %38*<<;<7<<7<=<<<;<<<<<
# seq1 279 C 23  A..T,,.,.,...,,,.,..... ;75&<<<<<<<<<=<<<9<<:<<
# test_ref = 'TTTAGAGCGC'

test_pileup = ['....,...,,,,.,...,',
               ',..,,,.c..C.,,,',
               'AAaAAaaaA.,AaaaAa',
               '.C.,.g,gG.G...,G.',
               '..,C,..C,,CCCC.C,c',
               '.,..tT,,.t',
               '',
               '.Taat..,.,...,',
               '.cc',
               'aA.,,A.,...,a.,Aa...,',
               ]

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(description='Naive variant detector using pileupreads notation')
        parser.add_argument('--referencefastafile', metavar='reference fasta file', type=str, nargs='+',
                            help='reference fasta file')
        parser.add_argument('--dnabasesfilename', metavar='dna bases file', type=str, nargs='+', help='valid dna bases')
        parser.add_argument('--pileupreads', metavar='reads pileupreads notation', type=str, nargs='+',
                            help='reads pileupreads notation')
        parser.add_argument('--pileupformat', metavar='pileupreads format', type=str, nargs='+',
                            help='pile up format notation')
        parser.add_argument('--minreaddepth', metavar='minimum read depth', type=int, nargs='+',
                            help='minimum read depth')
        parser.add_argument('--assemblymetadata', metavar='assembly information', type=str, nargs='+',
                            help='assembly information / chromosome lengths')
        parser.add_argument('--maxreadlength', metavar='max read lengths for md5 checksum generator', type=int,
                            nargs='+', help='maximum read lengths for md5 checksum generator')
        args = parser.parse_args()
        validityofassemblymetadatafile: bool
        validityofreferencefile: bool
        organism: str
        assembly: str
        chrLengthshash: dict[str, str]
        listofpileupreads: list[str]
        # Keep this line:
        (validityofassemblymetadatafile, organism, assembly, chrLengthshash) = chrlengthreader.readchromosomelengths(
            ''.join(args.assemblymetadata))
        if validityofassemblymetadatafile:
            dnabases: list[str]
            validityofreferencefile, reference, dnabases = refgenome.referencegenomeparser(''.join(args.referencefastafile),
                                                                                 ''.join(args.dnabasesfilename))

            if validityofreferencefile:
                # find maximum read length and generate md5 checksums till that length
                md5.randomMD5HashForDNAsequencegenerator(args.maxreadlength[0], "C:\\Users\\visu4\\PycharmProjects\\variantcallingfrompileup\\data\\md5checksum.txt")
                validityofpileupreadsfile: bool
                listofpileupreads: list[str]
                pileupnotationpattern: str
                validityofpileupreadsfile, listofpileupreads, pileupnotationpattern = pileup.pileupnotationreader(args.pileupformat, args.pileupreads)
                # if pileupreads are in right format, check for data integrity / samtools checker
                if validityofpileupreadsfile:
                    smtools.samtools_output_checker(listofpileupreads,
                                                    args.minreaddepth,
                                                    chrLengthshash, ''.join(dnabases),
                            pileupnotationpattern, "C:\\Users\\visu4\\PycharmProjects\\variantcallingfrompileup\\data\\md5checksum.txt")

        # TO DO: print the ref positions (1-10), and the variant call at each position
        # (Bonus: can you use list comprehension to print the results with one line of code?)
        # print(call_set)
        # find genes of significance

    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        # x, y = exception.args  # unpack args
        # print('x =', x)
        # print('y =', y)
