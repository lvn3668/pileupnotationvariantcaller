# Author: Lalitha Viswanathan
# Chromosome Length Reader (checks validity of file containing assembly information and chromosome lengths)
# Used to check validity of samtools mpileup output

import os
from typing import TextIO, Union


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
        raise Exception(filename, ' file is empty')
    else:
        print('File ', filename, ' is not empty')
        chrlengthsfilehandle: TextIO = open(''.join(filename), 'r')
        numberoflinesinfile: list[str] = chrlengthsfilehandle.readlines()
        # Strips the newline character
        for line in numberoflinesinfile:
            # Extract assembly information
            if line.startswith("assembly"):
                assembly: str = line.split(':')[1].strip()
            # Extract organism information
            elif line.startswith("organism"):
                organism: str = line.split(':')[1].strip()
            # Parse chromosome number and length
            elif len(line.strip()) > 1 and line.startswith("chr"):
                chromosomenumber: str
                chromosomelength: int
                (chromosomenumber, chromosomelength) = line.split(':')
                chromosomenumber = chromosomenumber.replace('chr', '')
                chromosomelengths[chromosomenumber.strip()] = str(chromosomelength).strip()

    if organism is not None and assembly is not None:
        return True, organism, assembly, chromosomelengths
    else:
        return False, None, None, None