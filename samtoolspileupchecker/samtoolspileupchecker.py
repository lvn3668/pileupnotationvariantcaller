def samtools_output_checker(pileupreads: list[str], min_depth: int, chrlengths: hash, dnabases: str,
                            pileupnotation: str) -> tuple[bool, Exception]:
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
        for string in pileupreads:
            # split as chromosome, 1-based coordinate, reference base, number of reads covering that position, reads themselves
            # min depth calculated for reads
            # QNAME, FLAG, RNAME, POS, MAPQ (displayed numerically), RNEXT, PNEXT.
            # Sanity check on samtools mpileup output file
            (chromosomenumber, positiononchromosome, referencebase, numberofreads, readstrings, qname, flag, rname, pos,
             cigar, mapq, rnext, pnext, tlen, SEQ, QUAL) = string.split('\s+')
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
                return False, Exception("Insufficient read depth")
            elif set(readstrings.split()).difference(pileupnotation.split()) > 0:
                return False, Exception("read basse strings contain invalid character")
            elif qname == None:
                return False, Exception("Query template Name absent; bam file in invalid format")
            elif flag == None:
                return False, Exception("Bitwise Flag absent")
            elif rname == None:
                return False, Exception("reference template name absent")
            elif pos == None:
                return False, Exception("1 base left most mapping position absent ; invalid BAM file")
            elif mapq == None:
                return False, Exception("Mapping Quality absent; Invalid BAM file")
            elif cigar == None:
                return False, Exception("CIGAR string absent")
            elif rnext == None:
                return False, Exception("reference for next read is absent")
            elif pnext == None:
                return False, Exception("Position of next read / mate pair is absent")
            elif tlen == None:
                return False, Exception("Observed template length is null")
            elif SEQ == None:
                return False, Exception("Observed sequence length is null")
            elif QUAL == None or QUAL < 33:
                return False, Exception("Sequence quality invalid")
            else:
                return True
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)
