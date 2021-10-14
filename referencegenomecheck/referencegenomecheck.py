import os
import re

# Author: Lalitha Viswanathan
# Project: Variant detection / Naive variant caller
from typing import TextIO


def referencegenomeparser(refgenomefilename: str, dnabasesfilename: str) -> tuple[bool, str, list[str]]:
    """
    :rtype: object
    :return:
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
            linesindnabasesfile: list[str] = dnabasesfilehandle.readlines()
            print(" Number of lines in ", dnabasesfilename, " is ", len(linesindnabasesfile))

            # Strips the newline character
            for line in linesindnabasesfile:
                if len(line.strip()) > 1:
                    linesindnabasesfile.remove(line)
                    raise Exception("DNA Bases file not in right format; One base per line")
            dnabasesfilehandle.close()
        ### Read dnabases file

        pattern = ''.join(linesindnabasesfile)

        ### Read reference genome file
        if os.stat(refgenomefilename).st_size == 0:
            print(refgenomefilename, ' File is empty')
            raise Exception(refgenomefilename, ' file is empty')
        else:
            print(refgenomefilename, ' File is not empty')
            refgenomefilehandle: TextIO = open(refgenomefilename, "r", encoding="utf-8")
            linesindnabasesfile: list[str] = refgenomefilehandle.readlines()

            for line in linesindnabasesfile:

                if re.search(r'['+pattern+']', line.strip()):
                    continue
                elif line.startswith('>'):
                    linesindnabasesfile.remove(line.strip())

            reference = ''.join(linesindnabasesfile).replace('\n', '')
            refgenomefilehandle.close()
            if len(reference) > 0 and len(linesindnabasesfile) > 0:
                print("Inside true ")
                return True, reference, linesindnabasesfile
            else:
                print("Inside false ")
                return False, reference, linesindnabasesfile

    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        # x = exception.args  # unpack args
        # print('x =', x)
