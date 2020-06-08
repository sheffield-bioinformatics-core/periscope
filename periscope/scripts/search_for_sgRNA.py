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
import snakemake

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

def open_bed(bed):
    """
    open bed file and return a bedtools object
    :param bed:
    :return:
    """
    bed_object = BedTool(bed)
    return bed_object

def set_tags(read,score,amplicon,read_class,orf):
    """
    assign tags to the given pysam object so that they can be added to an outpit bam file
    :param read: the pysam read
    :param score: the alignment score
    :param amplicon: the assigned amplicon
    :param read_class: the read class (i.e. sgRNA or gRNA)
    :param orf: the read orf
    :return: the pysam read object
    """
    read.set_tag('XS',score)
    read.set_tag('XA',amplicon)
    read.set_tag('XC',read_class)
    read.set_tag('XO', orf)

    return read

def setup_counts(primer_bed_object):
    # set up dictionary for normalisation
    # need to get all regions in bed and make into a dict
    # { 71: { total_reads: x, genomic_reads: y, sg_reads: {orf:z,orf2:k},'normalised_sgRNA': {orf:i,orf2:t} } }
    total_counts = {}
    for primer in primer_bed_object:
        amplicon = int(primer["Primer_ID"].split("_")[1])
        if amplicon not in total_counts:
            total_counts[amplicon] = {'pool': primer["PoolName"], 'total_reads': 0, 'genomic_reads': 0, 'sg_reads': {},
                                      'normalised_sgRNA': {}, 'normalised_gRNA': 0}
    return total_counts

def main(args):

    # read input bam file
    inbamfile = pysam.AlignmentFile(args.bam, "rb")
    # get bam header so that we can use it for writing later
    bam_header = inbamfile.header.copy().to_dict()
    # open output bam with the header we just got
    outbamfile = pysam.AlignmentFile(args.output_prefix + "_periscope.bam", "wb", header=bam_header)

    # get mapped reads
    print("hello")
    mapped_reads = pysam.idxstats(args.bam).split("\n")[0].split("\t")[2]

    dir=os.path.dirname(__file__)
    # open the orfs bed file
    orf_bed_object = open_bed(args.orf_bed)
    # open the amplicons bed file
    amplicon_bed_object= open_bed(args.amplicon_bed)
    # open the artic primer bed file
    primer_bed_object=read_bed_file(args.primer_bed)
    # set the output reads filename
    outfile_reads = args.output_prefix + "_periscope_reads.tsv"
    # set the output counts file name
    outfile_counts = args.output_prefix + "_periscope_counts.tsv"
    # set the output amplcion file name
    outfile_amplicons = args.output_prefix + "_periscope_amplicons.tsv"


    # add headers to these files
    file_reads = open(outfile_reads, "w")
    file_reads.write("sample\tread_id\tposition\tread_length\torf\tscore\tclass\tamplicon\n")

    total_counts = setup_counts(primer_bed_object)


    # get amplicon for each ORF
    orf_amplicons={}
    # we use pybedtools intersect here to intersect between the orf starts and the artic amplicons
    # in the case of E which actually overlaps two amplicons this overwirtes the entry for the 1st
    # E isn't far enough from the right of the 1st amplicon - 90bp and so we don't get sgRNAs from this amplicon
    # if amplcions hit pool two then they need the current amplicon and the next amplicon
    intersect = orf_bed_object.intersect(amplicon_bed_object,loj=True)
    for orf in intersect:
        print(orf)
        orf_name = orf[3]
        amplicon = int(orf[7].split("_")[1])
        if orf_name not in orf_amplicons:
            orf_amplicons[orf_name] = []


        if (amplicon % 2) == 0:
            # amplicon is even so you need the next amplicon too
            if amplicon not in orf_amplicons[orf_name]:
                orf_amplicons[orf_name].append(amplicon)
                if amplicon !=98:
                    orf_amplicons[orf_name].append(amplicon+1)
        else:
            if amplicon not in orf_amplicons[orf_name]:
                orf_amplicons[orf_name].append(amplicon)

    print(orf_amplicons)

    # for every read let's decide if it's sgRNA or not
    for read in inbamfile:

        # do some filtering
        # for some reason the odd read has no sequence
        if read.seq == None:
            print("%s read has no sequence" %
                  (read.query_name), file=sys.stderr)
            continue
        # filter chimeras
        # if len(read.seq) > args.max_read_length:
        #     print("%s skipped as too long" %
        #           (read.query_name), file=sys.stderr)
        #     continue
        # # filter short fragments like primers
        # if len(read.seq) < args.min_read_length:
        #     print("%s skipped as too short:" %
        #           (read.query_name), file=sys.stderr)
        #     continue
        # filter unmapped
        if read.is_unmapped:
            print("%s skipped as unmapped" %
                  (read.query_name), file=sys.stderr)
            continue
        # filter supplementary
        if read.is_supplementary:
            print("%s skipped as supplementary" %
                  (read.query_name), file=sys.stderr)
            continue

        # we are searching for the leader sequence
        search = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCT'
        # search for the sequence
        result = search_reads(read,search)
        # add orf location to result
        result["read_orf"] = check_start(orf_bed_object, read)

        # use artic code to find primers called "find_primers"
        # returns (1, 1, {'chrom': 'MN908947.3', 'start': 21357, 'end': 21386, 'Primer_ID': 'nCoV-2019_71_LEFT', 'PoolName': 'nCoV-2019_1', 'direction': '+'})
        # in the case of sgRNAs this code fails to find the correct primer for the + direction. This is actualy a clue that
        # this read is an sgRNA
        # get the left primer
        left_primer = find_primer(primer_bed_object,read.reference_start,'+')
        # get the left primer amplicon (we don't actually use this)
        left_amplicon = int(left_primer[2]['Primer_ID'].split("_")[1])
        # get the right primer
        right_primer = find_primer(primer_bed_object, read.reference_end, '-')
        # get the right primer amplicon
        right_amplicon = int(right_primer[2]['Primer_ID'].split("_")[1])

        # add to the total count of the right amplicon
        total_counts[right_amplicon]["total_reads"] += 1

        # if the read starts in an orf - let's just use all reads for now?
        # if result["read_orf"] != None:
            # if the alignment score is high
            # then it's a sgRNA
        if result["align_score"] > int(args.score_cutoff):
            # assign the read class
            read_class = "sgRNA"
            # # add one to the total count for the right amplicon for sgRNA reads only
            # if result["read_orf"] not in total_counts[right_amplicon]["sg_reads"]:
            #     total_counts[right_amplicon]["sg_reads"][result["read_orf"]] = 0
            #
            # attempt to capture novel:
            if result["read_orf"] is None:
                result["read_orf"] = "novel_" + str(read.pos)
            # if it has been an assigned an orf add to the total_counts
            if result["read_orf"] not in total_counts[right_amplicon]["sg_reads"]:
                total_counts[right_amplicon]["sg_reads"][result["read_orf"]] = 0
            total_counts[right_amplicon]["sg_reads"][result["read_orf"]] += 1


        # if not then it's gRNA
        else:
            read_class = "gRNA"
            # add one to the total count for the right amplicon for rRNA reads only
            total_counts[right_amplicon]["genomic_reads"] += 1

        # set the tags of the read for the bam file
        set_tags(read,result["align_score"],right_amplicon,read_class,result["read_orf"])
        # write the read to the bam file
        outbamfile.write(read)


        # output the read and some details to our reads file
        output = args.sample+"\t"+result["read_id"] + "\t" + str(result["read_position"]) +"\t" + str(len(read.seq)) + "\t" + str(result["read_orf"]) + "\t" + str(result["align_score"])+ "\t" +read_class+ "\t" + str(right_amplicon)
        # write the line to our file
        file_reads.write(output+"\n")

    # we're done with reads now
    file_reads.close()
    outbamfile.close()
    pysam.index(args.output_prefix + "_periscope.bam")

    # now do normalisation, for every orf get count
    # some sgRNAs hit pool 2 which doesn't have the leader primer
    # these have two sets of evidence, sgRNA reads resulting from priming within the amplicon from 3' with no 5' primer
    # the 2nd is priming in pool 1 from primer1 to the 3' end of the 3' most pool 1 amplicon
    # we have this data in total counts, as each amplicon has counts of sgRNA and counts of g, we just need to be sure to use both

    print(total_counts)
    result = {}


    for amplicon in total_counts:

        if total_counts[amplicon]["genomic_reads"] > 0:
            for orf in total_counts[amplicon]["sg_reads"]:


                if orf not in orf_amplicons:
                    orf_amplicons[orf]=[amplicon]

                total_counts[amplicon]["normalised_sgRNA"][orf] = total_counts[amplicon]["sg_reads"][orf] / total_counts[amplicon]["genomic_reads"]


                total_counts[amplicon]["normalised_gRNA"] = total_counts[amplicon]["genomic_reads"] / int(mapped_reads)

    file_amplicons = open(outfile_amplicons, "w")
    file_amplicons.write("sample\tamplicon\traw_gRNA\ttotal_reads\torf\traw_sgRNA\tnormalised_sgRNA\tnormalised_gRNA\n")
    for amplicon in total_counts:

        if len(total_counts[amplicon]["sg_reads"])>0:
            for orf in total_counts[amplicon]["sg_reads"]:
                line = []
                line.append(args.sample)
                line.append(str(amplicon))
                line.append(str(total_counts[amplicon]["genomic_reads"]))
                line.append(str(total_counts[amplicon]["total_reads"]))
                line.append(str(orf))
                line.append(str(total_counts[amplicon]["sg_reads"][orf]))
                line.append(str(total_counts[amplicon]["normalised_sgRNA"][orf]))

                line.append(str(total_counts[amplicon]["normalised_gRNA"]))
                file_amplicons.write("\t".join(line) + "\n")

        else:
            line = []
            line.append(args.sample)
            line.append(str(amplicon))
            line.append(str(total_counts[amplicon]["genomic_reads"]))
            line.append(str(total_counts[amplicon]["total_reads"]))
            line.append(str("NA"))
            line.append(str(0))
            line.append(str(0))
            line.append(str(total_counts[amplicon]["normalised_gRNA"]))
            file_amplicons.write("\t".join(line) + "\n")

    file_amplicons.close()



    # TODO: this is not working properly? Missing ORF8 and ORF10 because they don't get classified if there is no sgRNA above?
    # TODO: somtimes an ORF unexpecitidly is attributed to an amplicon
    for orf in orf_amplicons:
        print(orf)
        result[orf] = {
            "sg_reads":0,
            "genomic_reads":0,
            "normalised_sgRNA":0,
            "total_reads":0,
            "normalised_gRNA":0
        }
        for amplicon in orf_amplicons[orf]:
            print(amplicon)
            print(result[orf])
            if orf in total_counts[amplicon]["sg_reads"]:
                result[orf]["sg_reads"]+=total_counts[amplicon]["sg_reads"][orf]
                result[orf]["normalised_sgRNA"] += total_counts[amplicon]["normalised_sgRNA"][orf]
            else:
                result[orf]["sg_reads"]= 0
                result[orf]["normalised_sgRNA"] = 0
            result[orf]["genomic_reads"] += total_counts[amplicon]["genomic_reads"]
            result[orf]["total_reads"] += total_counts[amplicon]["total_reads"]

            result[orf]["normalised_gRNA"] += total_counts[amplicon]["normalised_gRNA"]




    # construct final counts file
    file_counts = open(outfile_counts, "w")
    file_counts.write("sample\torf\tnormalised_sgRNA\n")

    for orf in result:
        line=[]
        line.append(args.sample)
        line.append(str(orf))
        line.append(str(result[orf]["normalised_sgRNA"]))

        file_counts.write("\t".join(line)+"\n")

    # all done
    file_counts.close()
    return True



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




