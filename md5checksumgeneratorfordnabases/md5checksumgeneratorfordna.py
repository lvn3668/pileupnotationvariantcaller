import hashlib
from random import choice

# Author: Lalitha Viswanathan
# MD5 Checksum function
# Generates hash of md5 checksums for randomly generated DNA sequences from length 1 through specified length (100000)
# From SAMTOOLS SEQ header, the MD5 Value is validated as being present / absent in hash
# If present, valid samtools file else invalid samtools file\

# no dehashing to check if md5 dehashed is equivalent to the reference genome fragment for which the read information / alignment information
# is provided

def isValidMD5(string: str)-> bool:
    return string.matches("^[a-fA-F0-9]{32}$")

def weightedchoice(items: list[str]) -> object:  # this doesn't require the numbers to add up to 100
    """

    :param items:
    :return:
    """
    return choice("".join(x * y for x, y in items))


def generateRandomDNA(length: int) -> str:
    randomDNAstr_listofmd5hashes: list[str] = []
    DNA = ""
    for count in range(length):
        DNA += weightedchoice([("C", 10), ("G", 20), ("A", 40), ("T", 30)])
    #randomDNAstring[hashlib.md5(DNA)] = DNA
    randomDNAstr_listofmd5hashes.append(hashlib.md5(DNA))

    return randomDNAstr_listofmd5hashes

# Key is MD5 Hash and value is DNA String
def randomMD5HashForDNAsequencegenerator(maxlength: int) -> dict:
    if maxlength < 0 or maxlength < 100000:
        raise Exception ("Invalid max length for random DNA sequence generation")
    randomDNAsequencesWithMD5Hash: dict[str, str] = {}
    for length in range(1, maxlength):
        randomDNAsequencesWithMD5Hash.put(generateRandomDNA(length))

    return randomDNAsequencesWithMD5Hash

def checkIfValidMD5(md5: str, dictofmd5hashes: dict) -> bool:
    if md5 in dictofmd5hashes.keys():
        return True
    else:
        return False