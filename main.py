# This is a sample Python script.
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# from typing import List

import argparse
from typing import Any

import numpy as np

import md5checksumgeneratorfordnabases.md5checksumgeneratorfordnasequences as md5
import pileupfilereader.pileupnotationreader as pileup
import referencegenomecheck.referencegenomecheck as refgenome
import variantcaller.variantcaller as vc
import chromosomelengthreader.readassemblyinfo as chrlengthreader

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
# ACGTN	acgtn	Base is a mismatch to the reference base
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


def call_variants(pileupreads: list[str], min_depth: int) -> str:
    """
    :type pileupreads: object
    :param pileupreads:
    :param min_depth:
    :return:
    """
    # TODO: Make a variant call at each position on ref, using the pileupreads
    try:
        call_set: list[str] = []
        chromosomenumber: str = None
        positiononchromosome: int = None
        referencebase: str = None
        numberofreads: int = None
        readstrings: str = None
        qname: str = None
        flag: str = None
        rname: str = None
        pos: int = None
        cigar: str = None
        mapq: str = None
        rnext: str = None
        pnext: str = None
        tlen: str = None
        SEQ: str
        QUAL: str
        readlengths: dict[Any, Any] = {}
        readstrands = {}
        referenceskip = {}
        deletedbases = {}
        readquality = {}
        for stringofallreadsmappingtothatpositiononreference in pileupreads:
            chromosomenumber, positiononchromosome, referencebase, numberofreads, readstrings, qname, flag, rname, pos,
            cigar, mapq, rnext, pnext, tlen, SEQ, QUAL = stringofallreadsmappingtothatpositiononreference.split('\t')
            if len(readstrings) > min_depth:
                readnumber: int = 0
                reads = np.ones(len(readstrings))
                # readstrings
                # after caret, check if strand is given or ascii encoded mapping quality
                # if '+' isfollowed by number it is insertion of that many bases
                #
                while readnumber < len(readstrings):
                    # reference match on fwd strand (0)
                    if readstrings[readnumber] == ',':
                        reads[chromosomenumber][positiononchromosome][readnumber] = 0
                    # reference match on reverse strand (1)
                    elif readstrings[readnumber] == '.':
                        reads[chromosomenumber][positiononchromosome][readnumber] = 1
                    # Mismatch on FWD strand (2)
                    elif readstrings[readnumber] in ['A', 'T', 'G', 'C']:
                        reads[chromosomenumber][positiononchromosome][readnumber] = 2
                    elif readstrings[readnumber] in ['a', 't', 'g', 'c']:
                        reads[chromosomenumber][positiononchromosome][readnumber] = 3
                    # reference skip due to CIGAR string on forward strand
                    elif readstrings[readnumber] == '>':
                        referenceskip[chromosomenumber][positiononchromosome][readnumber] = 0
                    # reference skip due to CIGAR string on reverse strand
                    elif readstrings[readnumber] == '<':
                        referenceskip[chromosomelengthreader][positiononchromosome][readnumber] = 1
                    # read stops at this position
                    # 1 based indexing
                    # & , %33 (????)
                    elif readstrings[readnumber] == '$':
                        readlengths[chromosomenumber][positiononchromosome][readnumber] = positiononchromosome
                        # deleted bases on fwd strand (also on reverse unless flag specified in samtools call)
                    elif readstrings[readnumber] == '*':
                        deletedbases[chromosomenumber][positiononchromosome][readnumber] = 0
                    elif readstrings[readnumber] == '#':
                        deletedbases[chromosomenumber][positiononchromosome][readnumber] = 1
                    elif readstrings[readnumber] == '^':
                        readstrands[chromosomenumber][positiononchromosome][readnumber] = 'start'
                        readquality[chromosomenumber][positiononchromosome][readnumber] = ord(readstrings[readnumber+1])-33
                        # Move forward twice because start of read and its mapping quality has been parsed
                        readnumber = readnumber+1

            readnumber = readnumber + 1
            uniquebases, countsofuniquebases = np.unique(reads, return_counts=True)
            frequencies = np.asarray((uniquebases, countsofuniquebases)).T

            if vc.checkforreferenceallelecall(uniquebases):
                call_set.append("ref")
            elif vc.checkforaltallelecallinpileup(uniquebases):
                call_set.append("alt")
            elif vc.checkforalthomozygouscall(uniquebases, frequencies):
                call_set.append("[ACGT]hom")
            elif vc.checkifATGCHeterozygous(uniquebases, frequencies):
                call_set.append("[ACGT]het")
            elif vc.checkifATGClow(uniquebases, frequencies):
                call_set.append("[ACGT]low")
            else:
                call_set.append("nocall")
        return call_set
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


#def print_hi(name):
    """

#    :param name:
    """
    # Use a breakpoint in the code line below to debug your script.
 #   print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser(description='Naive variant detector using pileupreads notation')
        parser.add_argument('--referencefastafile', metavar='reference fasta file', type=str, nargs='+',
                        help='reference fasta file')
        parser.add_argument('--dnabasesfilename', metavar='dna bases file', type=str, nargs='+', help='valid dna bases')
        parser.add_argument('--pileupreads', metavar='reads pileupreads notation', type=str, nargs='+', help='reads pileupreads notation')
        parser.add_argument('--pileupformat', metavar='pileupreads format', type=str, nargs='+', help='pile up format notation')
        parser.add_argument('--minreaddepth', metavar='minimum read depth', type=int, nargs='+', help='minimum read depth')
        parser.add_argument('--assemblymetadata', metavar='assembly information', type=str, nargs='+',
                            help='assembly information / chromosome lengths')
        parser.add_argument('--maxreadlength', metavar='max read lengths for md5 checksum generator', type=int, nargs='+', help='maximum read lengths for md5 checksum generator')
        args = parser.parse_args()
        validityofassemblymetadatafile: bool = False
        validityofreferencefile: bool = False
        organism: str = None
        assembly: str = None
        chrLengthshash: dict[str, str] = None
        print(args.referencefastafile)
        print(args.assemblymetadata)
        # Keep this line:
        (validityofassemblymetadatafile, chrLengthshash, organism, assembly) = chrlengthreader.readchromosomelengths(''.join(args.assemblymetadata))
        if validityofassemblymetadatafile:
            print("Reference file ", ''.join(args.referencefastafile))
            validityofreferencefile, reference = refgenome.referencegenomeparser(''.join(args.referencefastafile), ''.join(args.dnabasesfilename))
            if validityofreferencefile:
                print("After reading reference genome")
                print("Inside md5 generator ", args.maxreadlength[0])
                md5checksumhash: dict[str, str] = md5.randomMD5HashForDNAsequencegenerator(args.maxreadlength[0])

            pileupreads: list[str] = pileup.pileupnotationreader(args.pileupformat, args.pileupreads)
            # find maximum read length and generate md5 checksums till that length

        #call_set = call_variants(test_pileup, reference, 5)

        # TO DO: print the ref positions (1-10), and the variant call at each position
        # (Bonus: can you use list comprehension to print the results with one line of code?)
        #print(call_set)
        # find genes of significance

    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        #x, y = exception.args  # unpack args
        #print('x =', x)
        #print('y =', y)