#!/usr/bin/env python3
from periscope import __version__
import argparse
import sys
import os
import snakemake


def main():

    parser = argparse.ArgumentParser(description='periscopre: Search for sgRNA reads in artic network SARS-CoV-2 sequencing data. A tool from Sheffield Bioinformatics Core/Florey Institute',usage='''periscope [options]''')
    parser.add_argument('--fastq-dir',dest='fastq_dir', help='the folder containing the raw pass demultiplexed fastqs from the artic protocol', default="resources/test.bam")
    parser.add_argument('--output-prefix',dest='output_prefix', help='Prefix of the output file',default="test")
    parser.add_argument('--score-cutoff',dest='score_cutoff', help='Cut-off for alignment score of leader (45)',default=45)
    parser.add_argument('--artic-primers', dest='artic_primers', help='artic network primer version used',
                        default="V1")
    parser.add_argument('--threads', dest='threads', help='number of threads',
                        default="1")
    parser.add_argument('-r', '--resources', dest='resources', help="the path to the periscope resources directory - whereever you cloned periscope into")
    parser.add_argument('-d', '--dry-run', action='store_true', help="perform a snakemake dryrun")
    parser.add_argument('-n', '--novel',dest='novel', action='store_true', help="don't restrict to known ORF sites",default=False)
    parser.add_argument('--max-read-length', dest='max_read_length', help="filter reads that have a length above this value (900)",default=900)
    parser.add_argument('--min-read-length', dest='min_read_length', help="filter reads that have a length below this value (200)",default=200)
    parser.add_argument('-f', '--force', action='store_true', help="Overwrite all output", dest="force")
    parser.add_argument('--sample', help='sample id', default="SHEF-D2BD9")

    args = parser.parse_args()

    # check if fastq_dir exists
    if not os.path.exists(args.fastq_dir):
        print("%s fastq directory must exist" % (args.fastq_dir), file=sys.stderr)
        exit(1)

    # check if version number is corect
    version = args.artic_primers.upper()
    if version not in ["V1", "V2", "V3"]:
        print("%s artic primer version incorrect" % (args.version), file=sys.stderr)
        exit(1)
    else:
        amplicons_bed="artic_amplicons_{}.bed".format(version)
        primers_bed="artic_primers_{}.bed".format(version)
        interest_bed = "artic_amplicons_of_interest.bed"


    # run snakemake pipeline 1st
    dir = os.path.join(os.path.dirname(__file__))
    scripts_dir= os.path.join(dir, 'scripts')




    config = {
        "fastq_dir": args.fastq_dir,
        "output_prefix": args.output_prefix,
        "scripts_dir":scripts_dir,
        "resources_dir":args.resources,
        "amplicon_bed": amplicons_bed,
        "interest_bed": interest_bed,
        "primer_bed": primers_bed,
        "orf_bed": "orf_start.bed",
        "score_cutoff": args.score_cutoff,
        "reference_fasta": 'nCoV-2019.reference.fasta',
        "sample": args.sample,
        "threads": args.threads,
        "min_read_length": args.min_read_length,
        "max_read_length": args.max_read_length
    }
    print(primers_bed)

    if args.novel == True:
        print("here")
        snakefile = os.path.join(scripts_dir, 'SnakefileN')
    else:
        print("this")
        snakefile = os.path.join(scripts_dir, 'Snakefile')
    print(snakefile)
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)
    else:
        print("Found the snakefile")

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=args.force, force_incomplete=True,
                                 config=config, cores=int(args.threads), lock=False
                                 )
    if status:  # translate "success" into shell exit code of 0
        exit(0)

    exit(1)






