import argparse
from .search_for_sgRNA import setup_counts







def reads_parser(read_file):

    with open(read_file, "r") as f:
        for line in f:
            sample,read_id,position,read_length,orf,score,read_class,amplicon=line.rstip().split("\t")
    # file looks like this
    # sample	read_id	position	read_length	orf	score	class	amplicon
    # SHEF-CAB89	acfa404a-736e-4c51-b809-648e0f5af7a4	14	497	None	46.0	gRNA	1



def main(args):
    reads_parser(args.reads_file)
    pass



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='periscopre: convert reads tsv into counts - useful for reprocessing with a different score')
    parser.add_argument('--output-prefix',dest='output_prefix', help='Prefix of the output file')
    parser.add_argument('--score-cutoff',dest='score_cutoff', help='Cut-off for alignment score of leader (45)',default=45)
    parser.add_argument('--sample', help='sample id')
    parser.add_argument('--reads-file',dest="reads_file", help='the reads file from your 1st periscope run')

    args = parser.parse_args()
    print(args.score_cutoff)

    main(args)