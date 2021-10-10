import re


def samtools_output_checker(pileupreads: list[str], min_depth: int, chrlengths: hash, dnabases: str,
                            pileupnotation: str) -> Union[tuple[bool, Exception], bool]:
    """
    :type pileupnotation: object
    :param pileupnotation: 
    :param dnabases: 
    :param chrlengths: 
    :type pileupreads: object
    :param pileupreads:
    :param min_depth:
    :return:
    """
    # TODO: Make a variant call at each position on ref, using the pileupreads
    try:

        chromosomenumber: str
        positionchromosome: int
        referencebase: str
        numberofreads: int
        readstrings: str
        qname: str
        flag: str
        rname: str
        pos: int
        cigar: str
        mapq: str
        rnext: str
        pnext: str
        tlen: str
        SEQ: str
        QUAL: str
        qnamepattern = re.compile('[!-?A-~]{1,254')
        cigarpattern = re.compile('\*|([0-9]+[MIDNSHIPX=1])+')
        seqpattern = re.compile('\*|[A-Za-z=.]+')
        qualpattern = re.match('[!-~]+')
        for string in pileupreads:
            # split as chromosome, 1-based coordinate, reference base, number of reads covering that position, reads themselves
            # min depth calculated for reads
            # QNAME, FLAG, RNAME, POS, MAPQ (displayed numerically), RNEXT, PNEXT.
            # Sanity check on samtools mpileup output file

            (chromosomenumber, positiononchromosome, referencebase, numberofreads, readstrings, qname, flag, rname, pos,
             cigar, mapq, rnext, pnext, tlen, SEQ, QUAL) = string.split('\t')
            if (chromosomenumber not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                                         '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X',
                                         'y']):
                return False, Exception("Chromosome information missing")
            # check if position is beyond length of chromosome
            elif (positiononchromosome == 0) or (positiononchromosome > chrlengths.get(chromosomenumber)):
                return False, Exception("Invalid position match on chromosome")
            # dna bases read earlier ; check if reference base matches one of the dna bases
            # else raise exception
            elif set(referencebase.split()).difference(dnabases.split()) > 0:
                return False, Exception("Invalid reference base")
            elif numberofreads < 0 or min_depth < 0 or numberofreads < min_depth:
                return False, Exception("Insufficient read depth; Minimum specified read depth ",min_depth)
            elif numberofreads != len(readstrings):
                return False, Exception("Number of reads does not match pileup reads ", numberofreads, " ", len(readstrings))
            elif set(readstrings.split()).difference(pileupnotation.split()) > 0:
                return False, Exception("read base strings contain invalid character")
            elif qname is None or not re.match(qnamepattern):
                return False, Exception("Query template Name absent; bam file in invalid format")
            elif flag is None or flag not in [0, 2 ^ 16 - 1]:
                return False, Exception("Bitwise Flag absent or not matching SAM Specifications")
            elif rname is None:
                return False, Exception("reference template name absent")
            elif pos is None or pos not in [0, 2 ^ 31 - 1]:
                return False, Exception("1 base left most mapping position absent or not matching SAM specifications; invalid BAM file")
            elif mapq is None or mapq not in [0, 2 ^ 8 - 1]:
                return False, Exception("Mapping Quality absent or matching SAM specifications; Invalid BAM file")
            elif cigar is None or not re.match(cigarpattern):
                return False, Exception("CIGAR string absent or not matching SAM specification for CIGAR string")
            elif rnext is None:
                return False, Exception("reference for next read is absent")
            elif pnext is None or pnext not in [0, 2 ^ 31 - 1]:
                return False, Exception("Position of next read / mate pair is absent or not as per SAM specifications")
            elif tlen is None or tlen not in [-2 ^ 31 + 1, 2 ^ 31 - 1]:
                return False, Exception("Observed template length is null or not as per SAM specifications")
            elif SEQ is None or not re.match(seqpattern):
                return False, Exception("Observed sequence length is null or not as per SAM specifications")
            elif QUAL is None or ord(QUAL)-33 < 33 or not re.match(qualpattern):
                return False, Exception("Sequence quality invalid or not as per SAM Specifications ")
            else:
                return True
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)
