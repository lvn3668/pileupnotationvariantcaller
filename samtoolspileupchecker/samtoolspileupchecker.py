import re
from typing import Pattern, Union
import md5checksumgeneratorfordnabases.md5checksumgeneratorfordnasequences as md5


# @HD	VN:1.0	SO:coordinate
# @SQ	SN:1	LN:249250621	AS:NCBI37	UR:file:/data/local/ref/GATK/human_g1k_v37.fasta	M5:1b22b98cdeb4a9304cb5d48026a85128
# @SQ	SN:2	LN:243199373	AS:NCBI37	UR:file:/data/local/ref/GATK/human_g1k_v37.fasta	M5:a0d9851da00400dec1098a9255ac712e
# @SQ	SN:3	LN:198022430	AS:NCBI37	UR:file:/data/local/ref/GATK/human_g1k_v37.fasta	M5:fdfd811849cc2fadebc929bb925902e5
# @RG	ID:UM0098:1	PL:ILLUMINA	PU:HWUSI-EAS1707-615LHAAXX-L001	LB:80	DT:2010-05-05T20:00:00-0400	SM:SD37743	CN:UMCORE
# @RG	ID:UM0098:2	PL:ILLUMINA	PU:HWUSI-EAS1707-615LHAAXX-L002	LB:80	DT:2010-05-05T20:00:00-0400	SM:SD37743	CN:UMCORE
# @PG	ID:bwa	VN:0.5.4
# @PG	ID:GATK TableRecalibration	VN:1.0.3471	CL:Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate, TileCovariate], default_read_group=null, default_platform=null, force_read_group=null, force_platform=null, solid_recal_mode=SET_Q_ZERO, window_size_nqs=5, homopolymer_nback=7, exception_if_no_tile=false, ignore_nocall_colorspace=false, pQ=5, maxQ=40, smoothing=1

# 1:497:R:-272+13M17D24M	113	1	497	37	37M	15	100338662	0	CGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG	0;==-==9;>>>>>=>>>>>>>>>>>=>>>>>>>>>>	XT:A:U	NM:i:0	SM:i:37	AM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37
# 19:20389:F:275+18M2D19M	99	1	17644	0	37M	=	17919	314	TATGACTGCTAATAATACCTACACATGTTAGAACCAT	>>>>>>>>>>>>>>>>>>>><<>>><<>>4::>>:<9	RG:Z:UM0098:1	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:4	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37
# 19:20389:F:275+18M2D19M	147	1	17919	0	18M2D19M	=	17644	-314	GTAGTACCAACTGTAAGTCCTTATCTTCATACTTTGT	;44999;499<8<8<<<8<<><<<<><7<;<<<>><<	XT:A:R	NM:i:2	SM:i:0	AM:i:0	X0:i:4	X1:i:0	XM:i:0	XO:i:1	XG:i:2	MD:Z:18^CA19
# 9:21597+10M2I25M:R:-209	83	1	21678	0	8M2I27M	=	21469	-244	CACCACATCACATATACCAAGCCTGGCTGTGTCTTCT	<;9<<5><<<<><<<>><<><>><9>><>>>9>>><>	XT:A:R	NM:i:2	SM:i:0	AM:i:0	X0:i:5	X1:i:0	XM:i:0	XO:i:1	XG:i:2	MD:Z:35

def parsepggrouplines(samtoolsdict: dict[str, str]):
    currentpgid = 0
    programgrouplines = line.split('\t')
    for entry in programgrouplines:
        (field, value) = entry.split(':')
        if field in ["ID", "PN", "CL", "PP", "DS", "VN"]:
            if field == "ID" and value in rgidentifierinsamtoolsheader.keys():
                return False, Exception(
                    "RG Id (Read Group Id) in @RG header in samtools file must be unique")
            elif field == "CL" and value is None:
                return False, Exception("Command Line field not conforming to SAMTools specification")
            elif field == "DS" and value is None:
                return False, Exception(
                    "DS  (Description) is not conforming to Samtools specifications")
            elif field == "VN" and value is None:
                return False, Exception(
                    "VN  (Version Number) is not conforming to Samtools specifications")
            elif field == "PN" and value is None:
                return False, Exception(
                    "PN Program Name not conforming to Samtools specifications")
            elif field != "ID":
                samtoolsdict['@PG'][currentpgid][field] = value
            # Valid Sequence Number field ; Hence add to dict
            elif field == "ID":
                pgidentifierinsamtoolsheader[value] = ""
                currentpgid = value
        else:
            return False, Exception('Invalid fields in Header field in SAM file')


def parsereadgrouplines(samtoolsdict: dict[str, str], readgrouplinesinsamtoolsheader: list[str],
                        rgidentifierinsamtoolsheader: list[str]):
    try:
        currentrgid = 0
        for line in readgrouplinesinsamtoolsheader:
            print(readgrouplinesinsamtoolsheader)
            readgrouplines: list[str] = line.split('\t')
            for entry in readgrouplines:
                print(entry)
                if ':' in entry:
                    (field, value) = entry.split(':')
                else:
                    continue
                if field in ["ID", "BC", "CN", "DS", "DS", "FO", "KS", "LB", "PG", "PI", "PL", "PM", "PU", "SM"]:
                    print(field,"\t",value)
                    if field == "ID" and value in rgidentifierinsamtoolsheader.keys():
                        return False, Exception(
                            "RG Id (Read Group Id) in @RG header in samtools file must be unique")
                    elif field == "BC" and value is None:
                        return False, Exception("Barcode in @RG Field not conforming to SAMTools specification")
                    elif field == "PG" and value is None:
                        return False, Exception(
                            "PG  (Program used for processing the read group) is not conforming to Samtools specifications")
                    elif field == "PI" and value is None:
                        return False, Exception(
                            "PI  (Predicted median insert size) is not conforming to Samtools specifications")
                    elif field == "CN" and value is None:
                        return False, Exception(
                            "Sequence center producing the reads is empty / not conforming to Samtools specifications")
                    # Not covering cases when PL is different from the below or is unknown (such sam files are marked as invalid; assumed sam file since 1-based indexing is checked for)
                    elif field == "PL" and value not in ['CAPILLARY', 'DNBSEQ', 'HELICOS', 'ILLUMINA', 'IONTORRENT',
                                                         'LS454', 'ONT', 'PACBIO', 'SOLID']:
                        return False, Exception(
                            "Platform / technology to produce the reads is empty / not conforming to Samtools specifications")
                    elif field == "PM" and value is None:
                        return False, Exception(
                            " Further details about the PM field in SAMTOOLS file not as per description")
                    elif field == "SM" and value is None:
                        return False, Exception(" SM / Sample field in SAMTOOLS file not as per description")
                    elif field == "DS" and value is None:
                        return False, Exception(" DS / Description field in SAMTOOLS file not as per description")
                    elif field == "DT" and not value.isDate():
                        return False, Exception(
                            "DT field (date when sequencing run was performed) is empty / not conforming to Samtools specifications")
                    elif field == "KS" and value is None:
                        return False, Exception(
                            "KS Field / array of nucleotide bases that correspond to the key sequence in each read in Sam file not as per specification")
                    elif field == "LB" and value is None:
                        return False, Exception("Library field in SAM file not as per specification")
                    elif field != "ID":
                        samtoolsdict['@RG'][currentrgid][field] = value
                    # Valid Sequence Number field ; Hence add to dict
                    elif field == "ID":
                        rgidentifierinsamtoolsheader[value] = ""
                        currentrgid = value
                else:
                    return False, Exception('Invalid fields in Header field in SAM file')
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


def parsesequencelines(samtoolsdict: dict[str, str], sequencelinesinsamtoolsheader: list[str],
                       sequencenamesinsamtoolsheader: dict[str, str], md5checksumdictfilename: str) -> Union[
    tuple[bool, Exception], tuple[bool, dict[str, str]]]:
    currentseqid = 0

    for entry in sequencelinesinsamtoolsheader:
        print(entry, "\t", " sequencelines ")
        fields = entry.strip().split("\t")
        for val in fields:
            print(val, " fields ")
            if ':' in val and not val.startswith("UR"):
                (field, value) = val.split(':')
                print("Field ", field, " VALUE ", value)
                if field in ["SN", "LN", "AH", "AS", "AN", "DS", "M5", "SP", "TP"]:
                    if field == "SN" and value in sequencenamesinsamtoolsheader.keys():
                        return False, Exception(
                            "SN Value (Sequence Id) in @SEQ header in samtools file must be unique")
                    elif field == "LN" and (int(value) < 0 or int(value) > 2 ** 31 - 1):
                        print("Inside false max LN Value ", 2 ** 31 - 1, " actual length ", value, " comparison ",
                              int(value) < 2 ** 31 - 1)
                        return False, Exception(
                            "LN  (Reference sequence length) is not conforming to Samtools specifications")
                    elif field == "LN" and (0 < int(value) < 2 ** 31 - 1):
                        print("Inside LN assignment ")
                    elif field == "AS" and value is None:
                        return False, Exception(
                            "AS / Assembly field is empty / not conforming to Samtools specifications")
                    # elif field == "UR" and not value.startswith("http"):
                    #    return False, Exception(
                    #    "UR field is empty / not conforming to Samtools specifications")
                    elif field == "SP" and value is None:
                        return False, Exception("SP field / species field in Sam file not as per specification")
                    elif field == "TP" and value not in ['linear', 'circular']:
                        return False, Exception("Topology field in SAM file not as per specification")
                    elif field == "DS" and value is None:
                        return False, Exception(" DS / Description field in SAMTOOLS file not as per description")
                    elif field == "M5" and (
                            md5.isValidMD5(value) is False or md5.checkIfValidMD5(value,
                                                                                  md5checksumdictfilename) is False):
                        print("Inside false ")
                        return False, Exception(
                            " MD5 Checksum invalid ")
                    elif field in ["SN", "LN", "AS", "M5"]:
                        print("Inside assignment clause ")
                        # sequencenamesinsamtoolsheader[value] = ""
                        # currentseqid = value
                else:
                    return False, Exception('Invalid fields in Header field in SAM file')
            else:
                print("Found UR field")
                continue
        print("before returning from function ")


def parseheaderlines(samtoolsdict: dict[str, str], headerlines: list[str],
                     formatinheaderlineversionpattern: Pattern[str],
                     subsortingalignmentspattern: Pattern[str]) -> Union[
    tuple[bool, Exception], tuple[bool, dict[str, str]]]:
    """

    :type samtoolsdict: object
    :param samtoolsdict:
    :param headerlines:
    :param formatinheaderlineversionpattern:
    :param subsortingalignmentspattern:
    :return:
    """
    print(headerlines.strip(), " Inside parse header lines function ")
    headerlinesarr = headerlines.strip().split("\t")
    headerdict = {"VN": [], "SO": [], "SS": [], "GO": []}

    print(headerlines)
    for entry in headerlinesarr:
        print(entry)
        if ':' in entry:
            (field, value) = entry.split(':')
        else:
            continue
        if field in ["VN", "SO", "SS", "GO"]:
            print(field, "\t", value, "\t", formatinheaderlineversionpattern)
            if (
                    field == "VN" and formatinheaderlineversionpattern.match(value)
                    or field == "SO" and value in ['unknown', 'unsorted', 'queryname', 'coordinate']
                    or field == "GO" and value in ['none', 'query', 'reference']
                    or field == "SS" and subsortingalignmentspattern.match(value)
            ):
                print("field ", field, " value ", value)
                headerdict[field].append(value)
                print("After assigning value ")
        else:
            return False, Exception(" Invalid header field in samtools: ", field, " value ", value)

    samtoolsdict['@HD'].append(headerdict)
    print("Before returning samtoolsdict")
    return True, samtoolsdict


def parseseqalignlines(samtoolsdict: dict[str, str], sequencelinesinsamtoolsheader: list[str],
                       sequencenamesinsamtoolsheader: dict[str, str], md5checksumdict: dict[str, str],
                       chrlengths: dict[str, str], min_depth: int, dnabases: str,
                       pileupnotation: str, qnamepattern: Pattern, cigarpattern: Pattern, seqpattern: Pattern,
                       qualpattern: Pattern):
    for line in sequencelinesinsamtoolsheader:
        print(line)
        (chromosomenumber, positiononchromosome, referencebase, numberofreads, readstrings, qname, flag, rname, pos,
         cigar, mapq, rnext, pnext, tlen, SEQ, QUAL) = line.split('\t')
        if (chromosomenumber not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                                     '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X',
                                     'y']):
            print("Inside chromosome number check ")
            return False, Exception("Chromosome information missing")
        # check if position is beyond length of chromosome
        elif (positiononchromosome == 0) or (positiononchromosome > chrlengths.get(chromosomenumber)):
            return False, Exception("Invalid position match on chromosome")
        # dna bases read earlier ; check if reference base matches one of the dna bases
        # else raise exception
        elif set(referencebase.split()).difference(dnabases.split()) > 0:
            return False, Exception("Invalid reference base")
        elif numberofreads < 0 or min_depth < 0 or numberofreads < min_depth:
            return False, Exception("Insufficient read depth; Minimum specified read depth ", min_depth)
        elif numberofreads != len(readstrings):
            return False, Exception("Number of reads does not match pileup reads ", numberofreads, " ",
                                    len(readstrings))
        elif set(readstrings.split()).difference(pileupnotation.split()) > 0:
            return False, Exception("read base strings contain invalid character")
        elif qname is None or not qnamepattern.match(qname):
            return False, Exception("Query template Name absent; bam file in invalid format")
        elif flag is None or flag not in [0, 2 ^ 16 - 1]:
            return False, Exception("Bitwise Flag absent or not matching SAM Specifications")
        elif rname is None:
            return False, Exception("reference template name absent")
        elif pos is None or pos not in [0, 2 ^ 31 - 1]:
            return False, Exception(
                "1 base left most mapping position absent or not matching SAM specifications; invalid BAM file")
        elif mapq is None or mapq not in [0, 2 ^ 8 - 1]:
            return False, Exception("Mapping Quality absent or matching SAM specifications; Invalid BAM file")
        elif cigar is None or not cigarpattern.match(cigar):
            return False, Exception("CIGAR string absent or not matching SAM specification for CIGAR string")
        elif rnext is None:
            return False, Exception("reference for next read is absent")
        elif pnext is None or pnext not in [0, 2 ^ 31 - 1]:
            return False, Exception("Position of next read / mate pair is absent or not as per SAM specifications")
        elif tlen is None or tlen not in [-2 ^ 31 + 1, 2 ^ 31 - 1]:
            return False, Exception("Observed template length is null or not as per SAM specifications")
        elif SEQ is None or not seqpattern.match(SEQ):
            return False, Exception("Observed sequence length is null or not as per SAM specifications")
        elif QUAL is None or ord(QUAL) - 33 < 33 or not qualpattern.match(QUAL) or ord(QUAL) - 33 == 255:
            return False, Exception("Sequence quality invalid or not as per SAM Specifications ")
        else:
            return True


def samtools_output_checker(pileupreads: list[str], min_depth: int, chrlengths: hash, dnabases: str,
                            pileupnotation: str, md5checksumdict: dict) -> Union[tuple[bool, Exception], bool]:
    """
    :param md5checksumdict:
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
        filelevelmetadata: bool = False
        headerseqreadgroupprogramlinepattern: Pattern[str] = re.compile(
            r"^@(HD|SQ|RG|PG)(\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$")

        qnamepattern: Pattern[str] = re.compile(r"[!-?A-~]{1,254}")
        cigarpattern: Pattern[str] = re.compile(r"\*|([0-9]+[MIDNSHIPX=1])+")
        seqpattern: Pattern[str] = re.compile(r"\*|[A-Za-z=.]+")
        qualpattern: Pattern[str] = re.compile(r"[!-~]+")
        samtoolsdict: dict[str, str] = {"@HD": [], "@SQ": [], "@RG": [], "@PG": []}

        rgidentifierinsamtoolsheader: dict[str, str] = {}
        pgidentifierinsamtoolsheader: dict[str, str] = {}
        validityofsamtoolsfile: bool
        sequencelinesinsamtoolsheader: list[str] = []
        readgrouplinesinsamtoolsheader: list[str] = []
        pggrouplinesinsamtoolsheader: list[str] = []
        commentgrouplinesinsamtoolsheader: list[str] = []
        sequencelinesinsamtoolsfile: list[str] = []
        for line in pileupreads:
            # split as chromosome, 1-based coordinate, reference base, number of reads covering that position, reads themselves
            # min depth calculated for reads
            # QNAME, FLAG, RNAME, POS, MAPQ (displayed numerically), RNEXT, PNEXT.
            # Sanity check on samtools mpileup output file
            print(line)
            # There can be only one HD Line which must conform to specified regex
            if line.startswith('@HD') and (
                    not re.match(headerseqreadgroupprogramlinepattern, line) or filelevelmetadata is True):
                return False, Exception("Invalid SAM File; Header field repeated multiple times")
            # if there is a HD line, and it matches the specified regex then parse the sub fields as present
            elif line.startswith('@HD') and (
                    re.match(headerseqreadgroupprogramlinepattern, line) and filelevelmetadata is False):
                # In HD Line, VN can be only of format major.minor version
                formatinheaderlineversionpattern: Pattern[str] = re.compile(r"^[0-9]+\.[0-9]+$")
                # sub sorting pattern can only be of format (type of sorting ):()
                subsortingalignmentspattern: Pattern[str] = re.compile(
                    r"(coordinate|queryname|unsorted)(:[A-Za-z0-9_-]+)+")
                validityofsamtoolsfile, samtoolsdict = parseheaderlines(samtoolsdict, line,
                                                                        formatinheaderlineversionpattern,
                                                                        subsortingalignmentspattern)
                print("After parsing header lines ", validityofsamtoolsfile)
            elif line.startswith('@SQ') and re.match(headerseqreadgroupprogramlinepattern, line):
                print("Inside SEQ line ", line)
                sequencenamesinsamtoolsheader: dict[str, str] = {}
                sequencelinesinsamtoolsheader.append(line)
            elif line.startswith('@SQ') and not re.match(headerseqreadgroupprogramlinepattern, line):
                return False, Exception("Sequence fields[@SQ] in samtools file not conforming to SAM specification")
            elif line.startswith('@RG') and re.match(headerseqreadgroupprogramlinepattern, line):
                readgrouplinesinsamtoolsheader.append(line)
            elif line.startswith('@RG') and not re.match(headerseqreadgroupprogramlinepattern, line):
                return False, Exception("Sequence fields[@RG] in samtools file not conforming to SAM specification")
            elif line.startswith('@PG'):
                pggrouplinesinsamtoolsheader.append(line)
            elif line.startswith('@CO'):
                commentgrouplinesinsamtoolsheader.append(line)
                commentgrouplines = line.split(
                    '\t')
                for entry in commentgrouplines:
                    (field, value) = entry.split(':')
                    samtoolsdict['@CO'][field] = value
            else:  # splitting sequence lines
                print("before splitting line ", line)
                sequencelinesinsamtoolsfile.append(line)

        print("Before parsing sequence lines ")
        if len(sequencelinesinsamtoolsheader) > 0:
            parsesequencelines(samtoolsdict, sequencelinesinsamtoolsheader,
                               sequencenamesinsamtoolsheader, md5checksumdict)
            print("After parsing sequence lines")
        if len(readgrouplinesinsamtoolsheader) > 0:
            print("Before read group lines in samtools header ")
            parsereadgrouplines(samtoolsdict, readgrouplinesinsamtoolsheader, rgidentifierinsamtoolsheader)
            print("After parsing read group lines")
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)
