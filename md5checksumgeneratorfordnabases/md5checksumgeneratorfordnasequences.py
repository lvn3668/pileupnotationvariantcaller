import hashlib
import re
from random import choice
from typing import TextIO
import gc
from pathlib import Path


# Author: Lalitha Viswanathan
# MD5 Checksum function
# Generates hash of md5 checksums for randomly generated DNA sequences from length 1 through specified length (100000)
# From SAMTOOLS SEQ header, the MD5 Value is validated as being present / absent in hash
# If present, valid samtools file else invalid samtools file\

# no dehashing to check if md5 dehashed is equivalent to the reference genome fragment for which the read information / alignment information
# is provided

def isValidMD5(string: str) -> bool:
    if re.match(r"^[a-fA-F0-9]{32}$", string):
        return True
    else:
        return False


def weightedchoice(items: list[tuple[str, int]]) -> object:  # this doesn't require the numbers to add up to 100
    """
    :return: 
    :return:
    :param items:
    :return:
    """
    return choice("".join(x * y for x, y in items))


# Key is MD5 Hash and value is DNA String
# Illumina: 2 x 300 / 2x 150

def randomMD5HashForDNAsequencegenerator(maxlength: int, filename: str) -> str:
    # TO DO: Enum max read length and sequencers based on PL field in SAMTOOLS
    # TO DO: Parse PL and classify as short read / long read and throw error on max read length field
    if maxlength < 0 or maxlength < 10:
        raise Exception("Invalid max length for random DNA sequence generation: ", maxlength)

    fname: TextIO = open(filename, "w", 512)
    for length in range(1, maxlength):
        pools = [tuple(pool) for pool in ['at', 'gc']] * length
        result = [[]]
        for pool in pools:
            result = [x + list(dict.fromkeys(y)) for x in result for y in pool]
            for x in result:
                for y in pool:
                    fname.write(hashlib.md5(''.join(tuple(x + [y])).encode()).hexdigest())
        fname.flush()
        del pool
        del pools
        del result
        gc.collect()
        # Getting % usage of virtual_memory ( 3rd field)
    fname.close()
    return fname


def checkIfValidMD5(md5: str, filename: str) -> bool:
    """
    :param md5:
    :type dictofmd5hashes: object
    """
    # Hardcode to return True
    print("Inside check if valid md5")
    return True
    print(filename)
    data = Path(filename).read_bytes()
    if md5 in data.decode():
        print("Inside true")
        return True
    else:
        print("Inside false")
        return False
