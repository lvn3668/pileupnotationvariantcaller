import os

# sequence id , 1 based coordinate, reference base, number of reads covering that base,
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


# Author: Lalitha Viswanathan
# Pileup read notation variant caller
from typing import TextIO, Union


def pileupnotationreader(pileupnotationfilename: str, pileupreadsfilename: str) -> Union[
    bool, tuple[bool, list[str], str]]:
    """
    :param pileupreadsfilename:
    :param pileupnotationfilename:
    """
    try:
        print("Pileup notation filename ", pileupnotationfilename)
        print("Pileup reads filename ", pileupreadsfilename)
        # check if size of file is 0
        if os.stat(''.join(pileupnotationfilename)).st_size == 0:
            print(pileupnotationfilename, ' File is empty')
            return False
        else:
            pileupnotationfilehandle: TextIO = open(''.join(pileupnotationfilename), 'r')
            pileuplines: list[str] = pileupnotationfilehandle.readlines()

            for line in pileuplines:
                if len(line) > 2:
                    pileuplines.remove(line)
                    return False
                else:
                    pileuplines[pileuplines.index(line)] = line.strip()

            pileupnotationpattern: str = ''.join(pileuplines)
            if os.stat(''.join(pileupreadsfilename)).st_size == 0:
                print(pileupreadsfilename, ' file is empty')
                return False
            else:
                print(pileupreadsfilename, ' file is not empty')
                pileupreadsfilehandle: TextIO = open(''.join(pileupreadsfilename), 'r')
                pileuplines: list[str] = pileupreadsfilehandle.readlines()
                for line in pileuplines:
                    if line is None:
                        pileuplines.remove(line)

                pileupreadsfilehandle.close()
        return True, pileuplines, pileupnotationpattern
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
