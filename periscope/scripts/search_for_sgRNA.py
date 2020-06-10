#!/usr/bin/env python3
from periscope import __version__

from Bio import pairwise2
import pysam
import argparse
from pybedtools import BedTool
import datetime
from artic.align_trim import find_primer
from artic.vcftagprimersites import read_bed_file
import sys
import os
import pprint as pp
import snakemake
from collections import namedtuple
class PeriscopeRead(object):
    def __init__(self, read):
        self.read = read



def get_mapped_reads(bam):
    mapped_reads = int(pysam.idxstats(bam).split("\n")[0].split("\t")[2])
    return mapped_reads

def check_start(bed_object,read):
    """
    find out where the read is in a bed file, in this case the ORF starts
    :param bed_object: bedtools object
    :param read: pysam read object
    :return: the orf
    """

    # reads with a pos of 0 make this fail so puting in try except works
    try:
        read_feature = BedTool(read.reference_name + "\t" + str(read.pos) + "\t" + str(read.pos), from_string=True)
        intersect = bed_object.intersect(read_feature)
        orf=intersect[0].name
        if len(intersect) > 1:
            print("odd")
    except:
        orf=None
    return orf

def search_reads(read,search):
    """
    given a pysam read object and a search string perform a localms alignment
    :param read: pysam read object
    :param search: DNA search string e.g. ATGTGCTTGATGC
    :return: dictionary containing the read_id, alignment score and the position of the read
    """
    align_score = pairwise2.align.localms(search, read.seq, 2, -2, -10, -.1,score_only=True)

    return {
        "read_id":  read.query_name,
        "align_score": align_score,
        "read_position": read.pos,
        "sequence": read.seq
    }


def find_amplicon(read,primer_bed_object):
    """
    use artic code to find primers called "find_primers"
    returns (1, 1, {'chrom': 'MN908947.3', 'start': 21357, 'end': 21386, 'Primer_ID': 'nCoV-2019_71_LEFT', 'PoolName': 'nCoV-2019_1', 'direction': '+'})
    in the case of sgRNAs this code fails to find the correct primer for the + direction. This is actualy a clue that
    this read is an sgRNA

    :param read:
    :param primer_bed_object:
    :return: the amplicon of the read
    """



    # get the left primer

    left_primer = find_primer(primer_bed_object, read.reference_start, '+')


    # get the right primer
    right_primer = find_primer(primer_bed_object, read.reference_end, '-')


    # get the left primer amplicon (we don't actually use this)
    left_amplicon = int(left_primer[2]['Primer_ID'].split("_")[1])
    # get the right primer amplicon
    right_amplicon = int(right_primer[2]['Primer_ID'].split("_")[1])

    # WARNING - LEFT_AMPLICON IS NOT RELIABLE FOR SG_RNA

    return dict(left_amplicon=left_amplicon,left_primer=left_primer,right_amplicon=right_amplicon,right_primer=right_primer)

def classify_read(read,score,score_cutoff,orf,amplicons):
    """
    classify read based on leader alignment score and other metrics
    :param score: the score
    :param score_cutoff: the user provided cut-off
    :return:
    """

    # some things I've learnt:
    # - if amplicons match it's more likely to by a gRNA but that doesn't hold true for reads that span amplicons - so score should still be 1st port of call
    # print(amplicons)
    quality=None
    if score > int(score_cutoff):
        quality = "HQ"
        if orf is not None:
            read_class = "sgRNA"
        else:
            read_class = "nsgRNA"

    elif score > 30:
        quality = "LQ"
        if orf is not None:
            read_class = "sgRNA"
        else:
            read_class = "nsgRNA"

    else:
        if orf is not None:
            quality = "LLQ"
            read_class = "sgRNA"
        else:
            read_class = "gRNA"

    # for those that have been classified as nsgRNA - do a final check - check not at amplicon edge
    # we see a lot of false positives at read ends

    if read_class == "nsgRNA":
        primer_start = amplicons["left_primer"][2]["start"]-5
        primer_end = amplicons["left_primer"][2]["end"]+5
        if primer_start <= read.pos <= primer_end:
            quality=None
            read_class="gRNA"

    if quality:
        print(read_class+"_"+quality)
        return read_class+"_"+quality
    else:
        print(read_class)
        return read_class


def open_bed(bed):
    """
    open bed file and return a bedtools object
    :param bed:
    :return:
    """
    bed_object = BedTool(bed)
    return bed_object


def setup_counts(primer_bed_object):
    # set up dictionary for normalisation
    # need to get all regions in bed and make into a dict
    # { 71: { total_reads: x, genomic_reads: y, sg_reads: {orf:z,orf2:k},'normalised_sgRNA': {orf:i,orf2:t} } }
    total_counts = {}
    for primer in primer_bed_object:
        amplicon = int(primer["Primer_ID"].split("_")[1])
        if amplicon not in total_counts:
            total_counts[amplicon] = {'pool': primer["PoolName"], 'total_reads': 0, 'gRNA': [], 'sgRNA_HQ': {}, 'sgRNA_LQ':{}, 'sgRNA_LLQ':{}, 'nsgRNA_HQ':{}, 'nsgRNA_LQ':{}}
    return total_counts


def calculate_normalised_counts(mapped_reads,total_counts,outfile_amplicon):
    with open(outfile_amplicon, "w") as f:
        header = ["sample", "amplicon", "mapped_reads", "orf", "quality", "gRNA_count", "gRPTH", "sgRNA_count", "sgRPHT",
              "sgRPTg"]
        f.write(",".join(header)+"\n")
        for amplicon in total_counts:

            # total count of gRNA for amplicon
            amplicon_gRNA_count = len(total_counts[amplicon]["gRNA"])

            # gRNA total count per 100l mapped reads
            amplicon_gRPTH = amplicon_gRNA_count / (mapped_reads / 100000)

            total_counts[amplicon]["gRPHT"] = {}

            for quality in ["HQ", "LQ", "LLQ"]:
                total_counts[amplicon]["sgRPHT_" + quality] = {}
                total_counts[amplicon]["sgRPTg_" + quality] = {}
                for orf in total_counts[amplicon]["sgRNA_" + quality]:
                    total_counts[amplicon]["gRPHT"][orf] = amplicon_gRPTH

                    amplicon_orf_sgRNA_count = len(total_counts[amplicon]["sgRNA_" + quality][orf])

                    # normalised per 100k total mapped reads
                    amplicon_orf_sgRPHT = amplicon_orf_sgRNA_count / (mapped_reads / 100000)

                    total_counts[amplicon]["sgRPHT_" + quality][orf] = amplicon_orf_sgRPHT

                    # normalised per 1000 gRNA reads from this amplicon
                    amplicon_orf_sgRPTg = amplicon_orf_sgRNA_count / (amplicon_gRNA_count / 1000)

                    total_counts[amplicon]["sgRPTg_" + quality][orf] = amplicon_orf_sgRPTg

                    line = []
                    line.append(args.sample)
                    line.append(str(amplicon))
                    line.append(str(mapped_reads))
                    line.append(str(orf))
                    line.append(str(quality))
                    line.append(str(amplicon_gRNA_count))
                    line.append(str(amplicon_gRPTH))
                    line.append(str(amplicon_orf_sgRNA_count))
                    line.append(str(amplicon_orf_sgRPHT))
                    line.append(str(amplicon_orf_sgRPTg))
                    f.write(",".join(line)+"\n")
    f.close()
    return total_counts


def summarised_counts_per_orf(total_counts,orf_bed_object):
    # TODO: would be nice if for ORFS even with 0 we get gRNA counts oututted
    # TODO: need to output novel counts
    # TODO: need to output novel counts
    result = {}
    for orf in orf_bed_object:


        if orf.name not in result:
            result[orf.name] = {}
            result[orf.name]["gRPHT"] = 0
            result[orf.name]["amplicons"] = []
            result[orf.name]["gRNA_count"] = 0
            for quality in ["LLQ", "LQ", "HQ"]:
                result[orf.name]["sgRNA_" + quality + "_count"] = 0

        for amplicon in total_counts:
            if orf.name in total_counts[amplicon]["gRPHT"]:
                result[orf.name]["gRPHT"] += total_counts[amplicon]["gRPHT"][orf.name]
                result[orf.name]["amplicons"].append(str(amplicon))
                result[orf.name]["gRNA_count"] += len(total_counts[amplicon]["gRNA"])

            for quality in ["LLQ", "LQ", "HQ"]:
                if orf.name in total_counts[amplicon]["sgRNA_" + quality]:
                    result[orf.name]["sgRNA_" + quality + "_count"] += len(total_counts[amplicon]["sgRNA_" + quality][orf.name])

                for metric in ["sgRPHT", "sgRPTg"]:
                    qmetric = metric + "_" + quality
                    if qmetric not in result[orf.name]:
                        result[orf.name][qmetric] = 0
                    if orf.name in total_counts[amplicon][qmetric]:
                        result[orf.name][qmetric] += total_counts[amplicon][qmetric][orf.name]
    return result

def output_summarised_counts(mapped_reads,result,outfile_counts):
    with open(outfile_counts,"w") as f:
        header = ["sample", "orf", "mapped_reads", "amplicons","gRNA_count", "sgRNA_HQ_count", "sgRNA_LQ_count", "sgRNA_LLQ_count", "gRHPT", "sgRPTg_HQ", "sgRPTg_LQ", "sgRPTg_LLQ",
                  "sgRPTg_ALL", "sgRPHT_HQ", "sgRPHT_LQ", "sgRPHT_LLQ", "sgRPHT_ALL"]
        f.write(",".join(header)+"\n")
        for orf in result:
            line = []
            line.append(args.sample)
            line.append(orf)
            line.append(str(mapped_reads))
            line.append("|".join(result[orf]["amplicons"]))
            line.append(str(result[orf]["gRNA_count"]))
            line.append(str(result[orf]["sgRNA_HQ_count"]))
            line.append(str(result[orf]["sgRNA_LQ_count"]))
            line.append(str(result[orf]["sgRNA_LLQ_count"]))

            line.append(str(result[orf]["gRPHT"]))
            line.append(str(result[orf]["sgRPTg_HQ"]))
            line.append(str(result[orf]["sgRPTg_LQ"]))
            line.append(str(result[orf]["sgRPTg_LLQ"]))
            sgRPTg_all = sum([result[orf]["sgRPTg_HQ"], result[orf]["sgRPTg_LQ"], result[orf]["sgRPTg_LLQ"]])
            line.append(str(sgRPTg_all))
            line.append(str(result[orf]["sgRPHT_HQ"]))
            line.append(str(result[orf]["sgRPHT_LQ"]))
            line.append(str(result[orf]["sgRPHT_LLQ"]))
            sgRPHT_all = sum([result[orf]["sgRPHT_HQ"], result[orf]["sgRPHT_LQ"], result[orf]["sgRPHT_LLQ"]])
            line.append(str(sgRPHT_all))

            f.write(",".join(line) + "\n")
        f.close()


def main(args):

    # read input bam file
    inbamfile = pysam.AlignmentFile(args.bam, "rb")
    # get bam header so that we can use it for writing later
    bam_header = inbamfile.header.copy().to_dict()
    # open output bam with the header we just got
    outbamfile = pysam.AlignmentFile(args.output_prefix + "_periscope.bam", "wb", header=bam_header)

    # get mapped reads
    mapped_reads = get_mapped_reads(args.bam)

    # open the orfs bed file
    orf_bed_object = open_bed(args.orf_bed)
    # open the artic primer bed file
    primer_bed_object=read_bed_file(args.primer_bed)
    # set the output reads filename
    outfile_reads = args.output_prefix + "_periscope_reads.tsv"
    # set the output counts file name

    # add headers to these files
    file_reads = open(outfile_reads, "w")
    file_reads.write("sample\tread_id\tposition\tread_length\torf\tscore\tclass\tamplicon\n")

    total_counts = setup_counts(primer_bed_object)
    count=0
    # for every read let's decide if it's sgRNA or not
    for read in inbamfile:
        count+=1
        print(count)
        if read.seq == None:
            print("%s read has no sequence" %
                  (read.query_name), file=sys.stderr)
            continue
        if read.is_unmapped:
            print("%s skipped as unmapped" %
                  (read.query_name), file=sys.stderr)
            continue
        if read.is_supplementary:
            print("%s skipped as supplementary" %
                  (read.query_name), file=sys.stderr)
            continue

        # find the amplicon for the read

        amplicons = find_amplicon(read, primer_bed_object)

        total_counts[amplicons["right_amplicon"]]["total_reads"] += 1


        # we are searching for the leader sequence
        search = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCT'

        # search for the sequence
        result = search_reads(read,search)

        # add orf location to result
        result["read_orf"] = check_start(orf_bed_object, read)

        # classify read based on prior information
        read_class = classify_read(read,result["align_score"],args.score_cutoff,result["read_orf"],amplicons)

        # store the attributes we have calculated with the read as tags
        read.set_tag('XS', result["align_score"])
        read.set_tag('XA', amplicons["right_amplicon"])
        read.set_tag('XC', read_class)
        read.set_tag('XO', result["read_orf"])


        # ok now add this info to a dictionary for later processing
        if "sgRNA" in read_class:
            if result["read_orf"] is None:
                result["read_orf"] = str(read.pos)

            if result["read_orf"] not in total_counts[amplicons["right_amplicon"]][read_class]:
                total_counts[amplicons["right_amplicon"]][read_class][result["read_orf"]] = []
            total_counts[amplicons["right_amplicon"]][read_class][result["read_orf"]].append(read)
        else:
            total_counts[amplicons["right_amplicon"]][read_class].append(read)

        # write the annotated read to a bam file
        outbamfile.write(read)

    outbamfile.close()
    pysam.index(args.output_prefix + "_periscope.bam")


    pp.pprint(total_counts)


    # go through each amplicon and do normalisations
    outfile_amplicons = args.output_prefix + "_periscope_amplicons.csv"
    total_counts = calculate_normalised_counts(mapped_reads,total_counts,outfile_amplicons)

    # summarise result into ORFs
    result = summarised_counts_per_orf(total_counts,orf_bed_object)

    # output summarised counts
    outfile_counts = args.output_prefix + "_periscope_counts.csv"
    output_summarised_counts(mapped_reads,result,outfile_counts)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='periscopre: Search for sgRNA reads in artic network SARS-CoV-2 sequencing data')
    parser.add_argument('--bam', help='bam file',default="resources/SHEF-D13D2_covid19-20200410-1586533428_all_raw_pass.bam")
    parser.add_argument('--output-prefix',dest='output_prefix', help='Prefix of the output file',default="test")
    parser.add_argument('--score-cutoff',dest='score_cutoff', help='Cut-off for alignment score of leader (45)',default=45)
    parser.add_argument('--orf-bed', dest='orf_bed', help='Cut-off for alignment score of leader (45)')
    parser.add_argument('--primer-bed', dest='primer_bed', help='Cut-off for alignment score of leader (45)')
    parser.add_argument('--amplicon-bed', dest='amplicon_bed', help='Cut-off for alignment score of leader (45)')
    parser.add_argument('--sample', help='sample id',default="SHEF-D2BD9")
    parser.add_argument('--max-read-length', dest='max_read_length', type=int, help="don't restrict to known ORF sites",
                        default=800)
    parser.add_argument('--min-read-length', dest='min_read_length', type=int,help="don't restrict to known ORF sites",
                        default=200)

    args = parser.parse_args()
    print(args.score_cutoff)

    periscope = main(args)
    if periscope:
        print("all done", file=sys.stderr)




