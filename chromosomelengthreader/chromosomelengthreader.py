# Author: Lalitha Viswanathan
# Chromosome Length Reader (checks validity of file containing assembly information and chromosome lengths)
# Used to check validity of samtools mpileup output

import os
from typing import TextIO


def readchromosomelengths(filename: str) -> tuple(bool, str, str, dict):
    """
    :return:
    :rtype: object
    :param filename:
    :return:
    """
    chromosomelengths: dict[str, str] = {}
    # check if size of file is 0
    if os.stat(''.join(filename)).st_size == 0:
        print('File is empty')
        raise Exception(filename, ' file is empty')
    else:
        print('File ', filename, ' is not empty')

        chrlengthsfilehandle: TextIO = open(''.join(filename), 'r')
        numberoflinesinfile: list[str] = chrlengthsfilehandle.readlines()
        print(" Number of lines ", len(numberoflinesinfile))

        # Strips the newline character
        for line in numberoflinesinfile:
            if line.startswith("assembly"):
                assembly = line.split(':')[1].strip()
                print(assembly)
            elif line.startswith("organism"):
                organism = line.split(':')[1].strip()
                print(organism)
            elif len(line.strip()) > 1 and line.startswith("chr"):
                chromosomenumber: str
                chromosomelength: int
                (chromosomenumber, chromosomelength) = line.split(':')
                chromosomelengths[chromosomenumber.strip()] = chromosomelength.strip()
            else:
                return False, Exception(" Chromosome lengths file not in right format")

    return True, organism, assembly, chromosomelengths
