import hashlib
import itertools
import zlib
from itertools import permutations
from random import choice


# Author: Lalitha Viswanathan
# MD5 Checksum function
# Generates hash of md5 checksums for randomly generated DNA sequences from length 1 through specified length (100000)
# From SAMTOOLS SEQ header, the MD5 Value is validated as being present / absent in hash
# If present, valid samtools file else invalid samtools file\

# no dehashing to check if md5 dehashed is equivalent to the reference genome fragment for which the read information / alignment information
# is provided
from typing import List


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

def randomMD5HashForDNAsequencegenerator(maxlength: int) -> dict:
    print("Inside random md5 hash generator")
    import io
    print("Default buffer size ", io.DEFAULT_BUFFER_SIZE)
    print("Inside random DNA sequence generator ")
    # TO DO: Enum max read length and sequencers based on PL field in SAMTOOLS
    # TO DO: Parse PL and classify as short read / long read and throw error on max read length field
    if maxlength < 0 or maxlength < 20:
        raise Exception("Invalid max length for random DNA sequence generation: ", maxlength)

    randomDNAsequencesWithMD5Hash: dict[str, str] = {}
    #fname = open(filename, "w", 512)
    maxlength  = 4
    for length in range(1, maxlength):
        print("Random DNA Sequence of ", length)

        # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
        print(length)
        pools = [tuple(pool) for pool in ['atgc', 'ATGC']] * length
        result = [[]]
        print("Pools ", pools)
        for pool in pools:
            print("Pool " , pool, " result ", result)
            result = [x + [y] for x in result for y in pool]
            print("Result ", result)

        for prod in result:
            op = yield tuple(prod)

        for item in op:
            print (item)
                #fname.write(str(hashlib.md5(item.encode()))+"\n")




        #for index, value in enumerate([''.join(x) for x in itertools.product('atgc', repeat=length)]):
            #randomDNAsequencesWithMD5Hash[hashlib.md5(value.encode())] = ""
            #fname.write(str(hashlib.md5(value.encode()))+"\n")
        #fname.flush()
    #fname.close()
    return randomDNAsequencesWithMD5Hash


def checkIfValidMD5(md5: str, dictofmd5hashes: dict) -> bool:
    """

    :param md5:
    :type dictofmd5hashes: object
    """
    if md5 in zlib.decompress(dictofmd5hashes.keys()):
        return True
    else:
        return False
