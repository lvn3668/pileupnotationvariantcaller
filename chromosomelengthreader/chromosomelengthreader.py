# Author: Lalitha Viswanathan
# Chromosome Length Reader (checks validity of file containing assembly information and chromosome lengths)
# Used to check validity of samtools mpileup output

import os
from types import Union
from typing import TextIO


def readchromosomelengths(filename: str) -> Union[tuple[bool, None, None, None], tuple[bool, str, str, dict[str, str]]]:
    """
    :return:
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
            # Extract assembly information
            if line.startswith("assembly"):
                assembly = line.split(':')[1].strip()
                print(assembly)
            # Extract organism information
            elif line.startswith("organism"):
                organism = line.split(':')[1].strip()
                print(organism)
            # Parse chromosome number and length
            elif len(line.strip()) > 1 and line.startswith("chr"):
                chromosomenumber: str
                chromosomelength: int
                (chromosomenumber, chromosomelength) = line.split(':')
                chromosomelengths[chromosomenumber.strip()] = str(chromosomelength).strip()
            else:
                return False, None, None, None

    return True, organism, assembly, chromosomelengths
