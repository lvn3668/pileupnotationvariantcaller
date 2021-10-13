import hashlib
import sys
import zlib
from random import choice
import io
from typing import TextIO

import psutil
import gc

# Author: Lalitha Viswanathan
# MD5 Checksum function
# Generates hash of md5 checksums for randomly generated DNA sequences from length 1 through specified length (100000)
# From SAMTOOLS SEQ header, the MD5 Value is validated as being present / absent in hash
# If present, valid samtools file else invalid samtools file\

# no dehashing to check if md5 dehashed is equivalent to the reference genome fragment for which the read information / alignment information
# is provided

def isValidMD5(string: str) -> bool:
    return string.matches("^[a-fA-F0-9]{32}$")


def weightedchoice(items: list[tuple[str, int]]) -> object:  # this doesn't require the numbers to add up to 100
    """
    :return: 
    :param items:
    :return:
    """
    return choice("".join(x * y for x, y in items))


# Key is MD5 Hash and value is DNA String
# Illumina: 2 x 300 / 2x 150

def randomMD5HashForDNAsequencegenerator(maxlength: int, filename: str):
    print('RAM memory % used:', psutil.virtual_memory()[2])
    print("Inside random md5 hash generator")
    print("Default buffer size ", io.DEFAULT_BUFFER_SIZE)
    # TO DO: Enum max read length and sequencers based on PL field in SAMTOOLS
    # TO DO: Parse PL and classify as short read / long read and throw error on max read length field
    if maxlength < 0 or maxlength < 20:
        raise Exception("Invalid max length for random DNA sequence generation: ", maxlength)

    fname: TextIO = open(filename, "w", 512)
    maxlength = 13
    for length in range(1, maxlength):
        print("Random DNA Sequence of length ", length)

        # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111

        pools = [tuple(pool) for pool in ['at', 'gc']] * length
        result = [[]]
        for pool in pools:
            result = [x + list(dict.fromkeys(y)) for x in result for y in pool]
            print("size of result set for nmer of length ", length, " is ", sys.getsizeof(result))
            for x in result:
                for y in pool:
                    fname.write(hashlib.md5(''.join(tuple(x + [y])).encode()).hexdigest())


        #for prod in result:
        #    fname.write(str(hashlib.md5((''.join(tuple(prod)).encode())))+"\n")

        fname.flush()
        del pool
        del pools
        del result
        #del prod
        gc.collect()
        # Getting % usage of virtual_memory ( 3rd field)
        print('RAM memory % used:', psutil.virtual_memory()[2])
        #if (psutil.virtual_memory()[2] > 52):
        #    exit(1);
        #    raise Exception("")
    fname.close()


def checkIfValidMD5(md5: str, dictofmd5hashes: dict) -> bool:
    """

    :param md5:
    :type dictofmd5hashes: object
    """
    if md5 in zlib.decompress(dictofmd5hashes.keys()):
        return True
    else:
        return False
