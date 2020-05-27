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
    parser.add_argument('-d', '--dry-run', action='store_true', help="perform a scnakemake dryrun")
    parser.add_argument('-n', '--novel', action='store_true', help="don't restrict to known ORF sites")
    parser.add_argument('-f', '--force', action='store_true', help="Overwrite all output", dest="force")
    parser.add_argument('--sample', help='sample id', default="SHEF-D2BD9")

    args = parser.parse_args()

    # check if fastq_dir exists

    # run snakemake pipeline 1st
    dir = os.path.join(os.path.dirname(__file__))
    scripts_dir= os.path.join(dir, 'scripts')
    resources_dir = os.path.join(dir, 'resources')

    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)
    else:
        print("Found the snakefile")

    version = args.artic_primers.upper()
    if version not in ["V1","V2","V3"]:
        exit(1)
    else:
        amplicons_bed="artic_amplicons_{}.bed".format(version)
        primers_bed="artic_primers_{}.bed".format(version)

    config = {
        "fastq_dir": args.fastq_dir,
        "output_prefix": args.output_prefix,
        "scripts_dir":scripts_dir,
        "resources_dir":args.resources,
        "amplicon_bed": amplicons_bed,
        "primer_bed": primers_bed,
        "orf_bed": "orf_start.bed",
        "score_cutoff": args.score_cutoff,
        "reference_fasta": 'nCoV-2019.reference.fasta',
        "sample": args.sample,
        "threads": args.threads
    }
    print(primers_bed)

    if args.novel:
        snakefile = os.path.join(scripts_dir, 'Snakefile')
    else:
        snakefile = os.path.join(scripts_dir, 'SnakefileN')

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=args.force, force_incomplete=True,
                                 config=config, cores=int(args.threads), lock=False
                                 )
    if status:  # translate "success" into shell exit code of 0
        exit(0)

    exit(1)






