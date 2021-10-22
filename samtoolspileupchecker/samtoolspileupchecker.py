import re
from typing import Pattern, Union
import md5checksumgeneratorfordnabases.md5checksumgeneratorfordnasequences as md5
import random


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
from main import call_variants
from variantcaller import variantcaller


def parsepggrouplines(samtoolsinfo: list[dict], pggrouplinesinsamtoolsheader: list[str],
                      pgidentifierinsamtoolsheader: dict[str, str]) -> object:
    """
    :param pggrouplinesinsamtoolsheader:
    :param pgidentifierinsamtoolsheader:
    """
    try:
        pggroupdict: dict[str, list] = {"ID": [], "PN": [], "CL": [], "PP": [], "DS": [], "VN": []}
        for entry in pggrouplinesinsamtoolsheader:
            pggrouplines: list[str] = entry.split('\t')
            for entry in pggrouplines:
                if ":" in entry:
                    (field, value) = entry.split(':')
                else:
                    continue
                if field in ["ID", "PN", "CL", "PP", "DS", "VN"]:
                    if (field == "ID" and value in pgidentifierinsamtoolsheader.keys()
                            or (field in ["CL", "DS", "VN", "PN", "ID"] and value is None)):
                        return False, Exception(
                            field, " in samtools file not as per samtools specification ")
                    elif field != "ID":
                        pggroupdict[field].append(value)
                    elif field == "ID":
                        pgidentifierinsamtoolsheader[value] = ""
                        pggroupdict[field].append(value)
                else:
                    return False

        samtoolsinfo.append(pggroupdict)
        return True, samtoolsinfo
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


def parsereadgrouplines(samtoolsinfo: list[dict], readgrouplinesinsamtoolsheader: list[str],
                        rgidentifierinsamtoolsheader: dict[str, str]):
    try:
        # @RG	ID:UM0098:1	PL:ILLUMINA	PU:HWUSI-EAS1707-615LHAAXX-L001	LB:80	DT:2010-05-05T20:00:00-0400	SM:SD37743	CN:UMCORE
        readgroupdict: dict[str, list] = {"ID": [], "BC": [], "CN": [], "DS": [], "F0": [], "KS": [], "LB": [],
                                          "PG": [],
                                          "PI": [], "PL": [], "PM": [], "PU": [], "SM": [], "DT": []}
        for line in readgrouplinesinsamtoolsheader:
            readgrouplines: list[str] = line.split('\t')
            for entry in readgrouplines:
                if ':' in entry and not entry.startswith("ID") and not entry.startswith("DT"):
                    (field, value) = entry.split(':')
                    if field in ["ID", "BC", "CN", "DS", "FO", "KS", "LB", "PG", "PI", "PL", "PM", "PU", "SM",
                                 "DT"] and value is not None:
                        readgroupdict[field].append(value)
                elif ':' in entry and entry.startswith("ID"):
                    (field, value, value2) = entry.split(':')
                    if value + "_" + value2 not in rgidentifierinsamtoolsheader.keys():
                        rgidentifierinsamtoolsheader[value + "_" + value2] = ""
                        readgroupdict[field].append(value)
                elif ':' in entry and entry.startswith("DT"):
                    m = re.compile("^DT:(\d{4}-\d{1,2}-\d{1,2})T(\d{1,2}:\d{1,2}:\d{1,2}-\d+)$")
                    g = m.search(entry)
                    if g:
                        field = "DT"
                        value = g.group(1)
                        value2 = g.group(2)
                        readgroupdict[field].append(value + "_" + value2)
        samtoolsinfo.append(readgroupdict)
        return True, samtoolsinfo
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)


def parsesequencelines(samtoolsinfo: list[dict], sequencelinesinsamtoolsheader: list[str],
                       sequencenamesinsamtoolsheader: dict[str, str], md5checksumdictfilename: dict[str, str]):
    """
    :param md5checksumdictfilename:
    :param sequencenamesinsamtoolsheader:
    :param samtoolsinfo:
    :rtype: object
    :type sequencelinesinsamtoolsheader: object
    """
    sequencedict: dict[str, list] = {"SN": [], "LN": [], "AH": [], "AS": [], "DS": [], "M5": [], "SP": [], "TP": [],
                                     "AN": []}
    for entry in sequencelinesinsamtoolsheader:
        fields = entry.strip().split("\t")
        for val in fields:
            if ':' in val and not val.startswith("UR"):
                (field, value) = val.split(':')
                if field in ["SN", "LN", "AH", "AS", "AN", "DS", "M5", "SP", "TP"]:
                    if field == "SN" and value in sequencenamesinsamtoolsheader.keys():
                        return False, Exception(
                            "SN Value (Sequence Id) in @SEQ header in samtools file must be unique")
                    elif field == "LN" and (int(value) < 0 or int(value) > 2 ** 31 - 1):
                        return False, Exception(
                            "LN  (Reference sequence length) is not conforming to Samtools specifications")
                    elif (
                            (field in ["AH", "AN", "M5", "AS", "SP", "TP", "DS"] and value is None)
                            or
                            (field == "MD5" and (
                                    md5.isValidMD5(value) is False or md5.checkIfValidMD5(value,
                                                                                          md5checksumdictfilename) is False))):
                        return False, Exception("Samtools field not as per specification ")
                    elif field == "SN":
                        sequencenamesinsamtoolsheader[value] = ""
                    elif field in ["AS", "AH", "AN", "DS", "M5", "SP", "TP"]:
                        sequencedict[field].append(value)
                else:
                    return False, Exception('Invalid fields in Header field in SAM file')
            else:
                continue
    samtoolsinfo.append(sequencedict)
    return True, samtoolsinfo


def parseheaderlines(samtoolsinfo: list[dict], headerlines: list[str],
                     formatinheaderlineversionpattern: Pattern[str],
                     subsortingalignmentspattern: Pattern[str]):
    """
    :param samtoolsinfo:
    :param headerlines:
    :param formatinheaderlineversionpattern:
    :param subsortingalignmentspattern:
    :return:
    """
    headerdict: dict[str, list] = {"VN": [], "SO": [], "SS": [], "GO": []}
    listofheadervals: list[str] = headerlines.strip().split("\t")
    for line in listofheadervals:
        if ':' in line:
            (field, value) = line.split(':')
        else:
            continue
        if field in ["VN", "SO", "SS", "GO"]:
            if (
                    field == "VN" and formatinheaderlineversionpattern.match(value)
                    or field == "SO" and value in ['unknown', 'unsorted', 'queryname', 'coordinate']
                    or field == "GO" and value in ['none', 'query', 'reference']
                    or field == "SS" and subsortingalignmentspattern.match(value)
            ):
                headerdict[field].append(value)
            else:
                return False

    samtoolsinfo.append(headerdict)
    return True, samtoolsinfo


def parseseqalignlines(samtoolsinfo: list[dict], sequencelinesinsamtoolsheader: list[str],
                       chrlengths: dict[str, str], min_depth: int, dnabases: str,
                       pileupnotation: str, qnamepattern: Pattern, cigarpattern: Pattern, seqpattern: Pattern,
                       qualpattern: Pattern, qname=None, flag=None, rname=None, cigar=None, tlen=None, mapq=None,
                       pnext=None, SEQ=None, QUAL=None, rnext=None, leftmostmappingpos=None):
    """

    :param samtoolsinfo: 
    :param sequencelinesinsamtoolsheader:
    :param chrlengths: 
    :param min_depth: 
    :param dnabases: 
    :param pileupnotation: 
    :param qnamepattern: 
    :param cigarpattern: 
    :param seqpattern: 
    :param qualpattern: 
    :return: 
    """
    seqaligndict: dict[str, list] = {"ALN": []}
    for line in sequencelinesinsamtoolsheader:
        (title, positiononchromosome, referencebase, numberofreads, readstrings, qual) = line.split('\t')
        # assign a random chromosome number
        # Assume human and not assigning Sex chromosomes
        randomnum: float = random.random()
        if randomnum < 0.5:
            chromosomenumber = random.randint(1, 23)
        elif randomnum < 0.7:
            chromosomenumber = 'X'
        else:
            chromosomenumber = 'Y'
        if chromosomenumber not in ['X', 'Y']:
            print(" chrlengths type ", type(chrlengths), " chromosome no ", chromosomenumber, " keys in dict ", chrlengths.keys(), " length of chromosome ", chrlengths.get(str(chromosomenumber)), " position on chromosome ", positiononchromosome)
        else:
            print(" chrlengths type ", type(chrlengths), " chromosome no ", chromosomenumber, " keys in dict ",
                  chrlengths.keys(), " length of chromosome ", chrlengths.get(chromosomenumber),
                  " position on chromosome ", positiononchromosome)
        print("*******chromosome number ", chromosomenumber, " position on chromosome ", positiononchromosome, " size of chromosome ", chrlengths.get(str(chromosomenumber)))
        print("reference base ", referencebase.lower().split(), " type of object ", type(referencebase.lower().split()))
        print("dna bases ", dnabases.lower().split())
        print(" Set ", set(''.join(referencebase.lower().split())) - set(''.join(dnabases.lower().split())))
        print(set(''.join(readstrings.split())))
        print(set(''.join(pileupnotation.split())))
        print("Difference between read strings and pile up notation ", set(''.join(readstrings.lower().split())) - set(''.join(pileupnotation.lower().split())) )
        if title is None:
            return False, None
        if chromosomenumber in [1, 2, 3, 4, 5, 6, 7,8,9,10, 11,12,13,14,15,16,17,18,19,20, 21,22,23,'X','Y']:
            print(" Valid chromosome number ")
            print(" type for position on chromosome ", positiononchromosome, " object type ", type(positiononchromosome))
            print(" chromosome length ", chrlengths.get(str(chromosomenumber)), "  type ", type(chrlengths.get(str(chromosomenumber))))
        # check if position is beyond length of chromosome
        if int(positiononchromosome) != 0:
            print(" position is valid")
        if positiononchromosome < str(chrlengths.get(str(chromosomenumber))):
            print(" position is valid ")
            # dna bases read earlier ; check if reference base matches one of the dna bases
            # else raise exception
            print(''.join(referencebase.lower().split()))
            print(''.join(dnabases.lower().split()))
            print("Set ", set(''.join(referencebase.lower().split())))
            print("Set ", set(''.join(dnabases.lower().split())))
        if len(set(''.join(referencebase.lower().split())) - set(''.join(dnabases.lower().split()))) == 0:
            print(" Valid check condition check for bases in reference base not in dna bases ")
            #return False, Exception("Invalid reference base")
        if int(numberofreads) > 0:
            print("clears num reads ; minimum depth ; required depth (cutoff to discard overlapping reads) ; stated number of reads and length of read string same ")
            print("minimum depth ", min_depth)
            print("type of object ", type(min_depth))
            print("joined as string ", int(min_depth[0]))
        if int(min_depth[0]) > 0:
            print("Valid check 2")
        if int(numberofreads) > int(min_depth[0]):
            print("Valid check 3 ", type(line), "\t", len(line), "\n")
            print("Number of reads ", int(numberofreads))
            print(" Length of line ", len(line))
        if int(numberofreads) == len(line):
            print("Valid check 4")
        #elif numberofreads != len(readstrings):
        #    return False, Exception("Number of reads does not match pileup reads ", numberofreads, " ",
        #                            len(readstrings))
        else:
            print("Inside else")
            seqaligndict["ALN"].append(line)
            print("adding to dict ")

        if qname is None or not qnamepattern.match(qname):
            print("Qname check is failing ")
            #return False, Exception("Query template Name absent; bam file in invalid format")
        if flag is None or flag not in [0, 2 ^ 16 - 1]:
            print("Flag check is failing ")
            #return False, Exception("Bitwise Flag absent or not matching SAM Specifications")
        if rname is None:
            print("Rname check is failing ")
            #return False, Exception("reference template name absent")
        if leftmostmappingpos is None or leftmostmappingpos not in [0, 2 ^ 31 - 1]:
            print(" left most mapping pos field is none ")
            #return False, Exception(
            #    "1 base left most mapping position absent or not matching SAM specifications; invalid BAM file")
        if mapq is None or mapq not in [0, 2 ^ 8 - 1]:
            print("mapq is null ")
            #return False, Exception("Mapping Quality absent or matching SAM specifications; Invalid BAM file")
        if cigar is None or not cigarpattern.match(cigar):
            print("cigar string is null ")
            #return False, Exception("CIGAR string absent or not matching SAM specification for CIGAR string")
        if rnext is None:
            print("rnext is null ")
            #return False, Exception("reference for next read is absent")
        if pnext is None or pnext not in [0, 2 ^ 31 - 1]:
            print("pnext is null")
            #return False, Exception("Position of next read / mate pair is absent or not as per SAM specifications")
        if tlen is None or tlen not in [-2 ^ 31 + 1, 2 ^ 31 - 1]:
            print("tlen is null ")
            #return False, Exception("Observed template length is null or not as per SAM specifications")
        if SEQ is None or not seqpattern.match(SEQ):
            print("seq is null")
            #return False, Exception("Observed sequence length is null or not as per SAM specifications")
        if QUAL is None or ord(QUAL) - 33 < 33 or not qualpattern.match(QUAL) or ord(QUAL) - 33 == 255:
            print("qual is null ")
            #return False, Exception("Sequence quality invalid or not as per SAM Specifications ")
        #else:
        seqaligndict["ALN"].append(line)
        print("adding to dict ****")

    samtoolsinfo.append(seqaligndict)
    print("Before returning samtoolsinfo ", len(seqaligndict.get("ALN")))

    return True, samtoolsinfo

def samtools_output_checker(pileupreads: list[str], min_depth: int, chrlengths: dict[str,str], dnabases: str,
                            pileupnotation: str, md5checksumdict: dict) -> Union[tuple[bool, Exception], bool]:
    """
    :type min_depth: object
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
        cigarpattern: Pattern[str] = re.compile(r"\*|([0-9]+[MIDNSHPX=1])+")
        seqpattern: Pattern[str] = re.compile(r"\*|[A-Za-z=.]+")
        qualpattern: Pattern[str] = re.compile(r"[!-~]+")
        samtoolsinfo: list[dict[str, str]] = []
        rgidentifierinsamtoolsheader: dict[str, str] = {}
        pgidentifierinsamtoolsheader: dict[str, str] = {}
        validityofsamtoolsfile: bool
        sequencelinesinsamtoolsheader: list[str] = []
        readgrouplinesinsamtoolsheader: list[str] = []
        pggrouplinesinsamtoolsheader: list[str] = []
        commentgrouplinesinsamtoolsheader: list[str] = []
        sequencelinesinsamtoolsfile: list[str] = []
        commentsdict: dict[str, str] = {}
        for line in pileupreads:
            # split as chromosome, 1-based coordinate, reference base, number of reads covering that position, reads themselves
            # min depth calculated for reads
            # QNAME, FLAG, RNAME, POS, MAPQ (displayed numerically), RNEXT, PNEXT.
            # Sanity check on samtools mpileup output file
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
                validityofsamtoolsfile, samtoolsinfo = parseheaderlines(samtoolsinfo, line,
                                                                        formatinheaderlineversionpattern,
                                                                        subsortingalignmentspattern)
            elif line.startswith('@SQ') and re.match(headerseqreadgroupprogramlinepattern, line):
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
                commentgrouplines: list[str] = line.split(
                    '\t')
                for entry in commentgrouplines:
                    if ':' in entry:
                        (field, value) = entry.split(':')
                        commentsdict[field].append(value)
            else:  # splitting sequence lines
                sequencelinesinsamtoolsfile.append(line)

        validityofsamtoolsfile: bool
        if len(sequencelinesinsamtoolsheader) > 0:
            validityofsamtoolsfile, samtoolsinfo = parsesequencelines(samtoolsinfo, sequencelinesinsamtoolsheader,
                                                                      sequencenamesinsamtoolsheader, md5checksumdict)
        if len(readgrouplinesinsamtoolsheader) > 0:
            validityofsamtoolsfile, samtoolsinfo = parsereadgrouplines(samtoolsinfo, readgrouplinesinsamtoolsheader,
                                                                       rgidentifierinsamtoolsheader)
            print("After parsing read group lines ", len(samtoolsinfo))
        if len(pggrouplinesinsamtoolsheader) > 0:
            validityofsamtoolsfile, samtoolsinfo = parsepggrouplines(samtoolsinfo, pggrouplinesinsamtoolsheader,
                                                                     pgidentifierinsamtoolsheader)
            print("After parsing pg group lines ", len(samtoolsinfo))
        if len(commentsdict) > 0:
            samtoolsinfo.append(commentsdict)
            print("After parsing comment dicts (no entries) ", len(samtoolsinfo))
        if len(sequencelinesinsamtoolsheader) > 0:
            validityofsamtoolsfile, samtoolsinfo = parseseqalignlines(samtoolsinfo, sequencelinesinsamtoolsfile,
                               chrlengths, min_depth, dnabases,
                               pileupnotation, qnamepattern, cigarpattern, seqpattern,
                               qualpattern)

            print("After parsing sequence alignment records ", len(samtoolsinfo))

        # Call variants
        call_variants(samtoolsinfo.get("ALN"), min_depth)
    except Exception as exception:
        print(type(exception))  # the exception instance
        print(exception.args)  # arguments stored in .args
        print(exception)  # __str__ allows args to be printed directly,
        x, y = exception.args  # unpack args
        print('x =', x)
        print('y =', y)
