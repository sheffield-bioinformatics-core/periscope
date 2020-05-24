from Bio import pairwise2
import pysam
import argparse
from pybedtools import BedTool
import datetime
from artic.align_trim import find_primer
from artic.vcftagprimersites import read_bed_file
import sys


def check_start(bed_object,read):
    """
    find out where the read is in a bed file, in this case the ORF starts
    :param bed_object: bedtools object
    :param read: pysam read object
    :return: the orf
    """
    read_feature = BedTool(read.reference_name+"\t"+str(read.pos)+"\t"+str(read.pos),from_string=True)
    intersect=bed_object.intersect(read_feature)
    try:
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

    # if align_score > 55:
    #     if read.is_reverse:
    #         print(">reverse\n"+read.seq)
    #     else:
    #         print(">forward\n" + read.seq)
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

def set_tags(read,score,amplicon,read_class):

    read.set_tag('XS',score)
    read.set_tag('XA',amplicon)
    read.set_tag('XC',read_class)
    read.set_tag('RG', read_class)

    return read


def main(args):

    inbamfile = pysam.AlignmentFile(args.bam, "rb")
    bam_header = inbamfile.header.copy().to_dict()
    outbamfile = pysam.AlignmentFile(args.output_prefix + "_periscope.bam", "wh", header=bam_header)

    orf_bed_object = open_bed('resources/orf_start.bed')
    amplicon_bed_object= open_bed('resources/artic_amplicons_V1.bed')
    primer_bed_object=read_bed_file('resources/artic_primers_V1.bed')
    outfile = args.output_prefix + "_seq_search.tsv"

    # set up dictionary for normalisation
    # need to get all regions in bed and make into a dict
    # { 71: { total_reads: x, genomic_reads: y, sg_reads: z } }
    total_counts= {}
    for primer in primer_bed_object:
        amplicon = int(primer["Primer_ID"].split("_")[1])
        if amplicon not in total_counts:
            total_counts[amplicon] = {'total_reads': 0, 'genomic_reads': 0, 'sg_reads': 0}

    # set up an expression summary - so per ORF?
    orf_counts={}
    for orf in orf_bed_object:
        print(orf.name)
        orf_counts[orf.name]=0

     # set up an expression summary - so per ORF?
    orf_norm = {}
    for orf in orf_bed_object:
        # print(orf.name)
        orf_norm[orf.name] = 0

    # get amplicon for each ORF
    orf_amplicons={}
    intersect = orf_bed_object.intersect(amplicon_bed_object,loj=True)
    for orf in intersect:
        orf_name = orf[3]
        amplicon = int(orf[7].split("_")[1])
        orf_amplicons[orf_name]=amplicon



    file = open(outfile, "w")
    file.write("sample\tread_id\tposition\torf\tscore\tsequence\n")
    print(datetime.datetime.now())




    for read in inbamfile:


        if read.is_unmapped:
            # print("%s skipped as unmapped" %
            #       (read.query_name), file=sys.stderr)
            continue
        if read.is_supplementary:
            # print("%s skipped as supplementary" %
            #       (read.query_name), file=sys.stderr)
            continue

        # print(read.query_name)
        search = 'AACCAACTTTCGATCTCTTGTAGATCTGTTCT'
        result = search_reads(read,search)
        # add orf location to result
        result["read_orf"] = check_start(orf_bed_object, read)

        # in the case of sgRNAs this code fails to find the correct primer for the + direction. This is actualy a clue that
        # this read is an sgRNA
        # get the left primer
        left_primer = find_primer(primer_bed_object,read.reference_start,'+')
        # get the left primer amplicon
        left_amplicon = int(left_primer[2]['Primer_ID'].split("_")[1])
        # get the right primer
        right_primer = find_primer(primer_bed_object, read.reference_end, '-')
        # get the right primer amplicon
        right_amplicon = int(right_primer[2]['Primer_ID'].split("_")[1])


        # (1, 1, {'chrom': 'MN908947.3', 'start': 21357, 'end': 21386, 'Primer_ID': 'nCoV-2019_71_LEFT', 'PoolName': 'nCoV-2019_1', 'direction': '+'})
        # if the left and right ampicon match then it's genome
        # todo: amplion 1 is genomic and contains leader, so here we need to do something else?
        if left_amplicon == right_amplicon:
            # print("%s ampilcon matches so genomic" %
            #       (read.query_name), file=sys.stderr)
            read_class = "gRNA"
            total_counts[right_amplicon]["total_reads"]+=1
            total_counts[right_amplicon]["genomic_reads"] += 1
        # if the left and right amplicon do not match it's sub-genomic
        elif left_amplicon != right_amplicon:
            # print("%s ampilcon mismatch so likely sub-genomic" %
            #       (read.query_name), file=sys.stderr)
            read_class = "sgRNA"
            total_counts[right_amplicon]["total_reads"] += 1
            total_counts[right_amplicon]["sg_reads"] += 1

        # ok so now, if read is an ORF, if leader is found, and artic is sgRNA then add to orf counts
        if result["read_orf"] != None:
            if result["align_score"] > int(args.score_cutoff):
                orf_counts[result["read_orf"]]+=1

        set_tags(read,result["align_score"],right_amplicon,read_class)

        # output = sample+"\t"+result["read_id"] + "\t" + str(result["read_position"]) + "\t" + str(result["read_orf"]) + "\t" + str(result["align_score"])+ "\t" + str(result["sequence"]) + "\n"
        output = args.sample+"\t"+result["read_id"] + "\t" + str(result["read_position"]) + "\t" + str(result["read_orf"]) + "\t" + str(result["align_score"])+ "\t" +read_class+ "\t" + str(right_amplicon)
        # print(output)
        outbamfile.write(read)
        file.write(output+"\n")


    print(total_counts)
    # now do normalisation, for every orf get count
    for orf in orf_counts:
        # get the amplicon to which the ORF belongs
        amplicon = orf_amplicons[orf]
        # get the total reads for that amplicon
        total_for_amplicon = total_counts[amplicon]["total_reads"]
        # exclude 0 counts
        if total_for_amplicon != 0:
            # divide orf sgRNA count by the total reads for that amplicon
            orf_norm[orf]=orf_counts[orf]/total_for_amplicon
        else:
            orf_norm[orf] = 'NA'
    print(orf_counts)
    print(orf_norm)
    print(datetime.datetime.now())
    file.close()
    outbamfile.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='periscopre: Search for sgRNA reads in artic network SARS-CoV-2 sequencing data')
    # parser.add_argument('--bam', help='bam file',default="resources/SHEF-C0E7B_covid19-20200421-1587480569_all_raw_pass.bam")
    parser.add_argument('--bam', help='bam file',default="/mnt/shared/covid19seq/Shared/playground/trs/test.sam")
    parser.add_argument('--output-prefix',dest='output_prefix', help='Prefix of the output file',default="test2")
    parser.add_argument('--score-cutoff',dest='score_cutoff', help='Cut-off for alignment score of leader (55)',default="55")

    parser.add_argument('--sample', help='sample id',default="SHEF-D2BD9")

    args = parser.parse_args()


    main(args)




