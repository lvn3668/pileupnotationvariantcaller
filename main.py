# This is a sample Python script.
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# from typing import List

import argparse
import numpy as np
import pileupfilereader.pileupnotationreader as pileup
import referencegenomecheck.referencegenomecheck as refgenome
import variantcaller.variantcaller as vc

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
from chromosomelengthreader import chromosomelengthreader

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
        chromosomenumber: str
        positiononchromosome: int
        referencebase: str
        numberofreads: int
        readstrings: str
        qname: str
        flag: str
        rname: str
        pos: int
        cigar: str
        mapq: str
        rnext: str
        pnext: str
        tlen: str
        SEQ: str
        QUAL: str
        readlengths = {}
        readstrands = {}
        referenceskip = {}
        deletedbases = {}
        for stringofallreadsmappingtothatpositiononreference in pileupreads:
            chromosomenumber, positiononchromosome, referencebase, numberofreads, readstrings, qname, flag, rname, pos,
            cigar, mapq, rnext, pnext, tlen, SEQ, QUAL = stringofallreadsmappingtothatpositiononreference.split('\s+')
            if len(stringofallreadsmappingtothatpositiononreference) > min_depth:
                readnumber = 0
                reads = np.ones(len(stringofallreadsmappingtothatpositiononreference))
                # readstrings
                # after caret, check if strand is given or ascii encoded mapping quality
                # if '+' isfollowed by number it is insertion of that many bases
                #
                while readnumber < len(stringofallreadsmappingtothatpositiononreference):
                    # reference match on fwd strand (0)
                    if stringofallreadsmappingtothatpositiononreference[readnumber] == ',':
                        reads[chromosomenumber][readnumber] = 0
                    # reference match on reverse strand (1)
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] == '.':
                        reads[chromosomenumber][readnumber] = 1
                    # Mismatch on FWD strand (2)
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] in ['A', 'T', 'G', 'C']:
                        reads[chromosomenumber][readnumber] = 2
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] in ['a', 't', 'g', 'c']:
                        reads[chromosomenumber][readnumber] = 3
                    # reference skip due to CIGAR string on forward strand
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] == '>':
                        referenceskip[chromosomenumber][readnumber][positiononchromosome] = 0
                    # reference skip due to CIGAR string on reverse strand
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] == '<':
                        referenceskip[chromosomelengthreader][readnumber][positiononchromosome] = 1
                    # read stops at this position
                    # 1 based indexing
                    # & , %33 (????)
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] == '$':
                        readlengths[chromosomenumber][readnumber][positiononchromosome] = 'end'
                        # deleted bases on fwd strand (also on reverse unless flag specified in samtools call)
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] == '*':
                        deletedbases[chromosomenumber][readnumber][positiononchromosome] = 0
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] == '#':
                        deletedbases[chromosomenumber][readnumber][positiononchromosome] = 1
                    elif stringofallreadsmappingtothatpositiononreference[readnumber] == '^':
                        readstrands[chromosomenumber][readnumber][positiononchromosome] = 'start'

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
        args = parser.parse_args()

        print(args.referencefastafile)
        print(args.assemblymetadata)
        # Keep this line:
        (Chrlengthshash, organism, assembly) = chromosomelengthreader.readchromosomelengths(args.assemblymetadata)
        reference = refgenome.referencegenomeparser(args.referencefastafile, args.dnabasesfilename)
        pileupreads: list[str] = pileup.pileupnotationreader(args.pileupformat, args.pileupreads)
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