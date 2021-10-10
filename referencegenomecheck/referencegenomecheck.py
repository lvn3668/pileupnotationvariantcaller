import os
import re

# Author: Lalitha Viswanathan
# Project: Variant detection / Naive variant caller
from typing import TextIO


def referencegenomeparser(refgenomefilename: str, dnabasesfilename: str) -> object:
    """
    :param dnabasesfilename:
    :type refgenomefilename: object
    """
    try:
        global dnabases
        dnabases = []
        reference: str = ""
        # check if size of file is 0
        if os.stat(''.join(dnabasesfilename)).st_size == 0:
            print('File is empty')
            raise Exception(dnabasesfilename, ' file is empty')
        else:
            print('File ', dnabasesfilename, ' is not empty')

            dnabasesfilehandle: TextIO = open(dnabasesfilename, 'r')
            numberoflinesinfile = dnabasesfilehandle.readlines()
            print(" Number of lines ", len(numberoflinesinfile))

            # Strips the newline character
            for line in numberoflinesinfile:
                if len(line.strip()) > 1:
                    numberoflinesinfile.remove(line)
                    raise Exception("DNA Bases file not in right format; One base per line")

        dnabasesfilehandle.close()
        if os.stat(refgenomefilename).st_size == 0:
            print(refgenomefilename, ' File is empty')
            raise Exception(refgenomefilename, ' file is empty')
        else:
            print(refgenomefilename, ' File is not empty')
            refgenomefilehandle: TextIO = open(refgenomefilename, "r", encoding="utf-8")
            lines: list[str] = refgenomefilehandle.readlines()

            for line in lines:
                if line.startswith('>') or re.search(r'[ATGCUatgcu]', line.strip()):
                    continue
                else:
                    lines.remove(line.strip())
            reference = ''.join(lines).replace('\n', '')
            refgenomefilehandle.close()

    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        # x = exception.args  # unpack args
        # print('x =', x)

    return reference
