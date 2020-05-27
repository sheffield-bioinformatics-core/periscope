import argparse








def reads_parser(read_file):
    pass


def main():
    pass



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='periscopre: convert reads tsv into counts - useful for reprocessing with a different score')
    parser.add_argument('--output-prefix',dest='output_prefix', help='Prefix of the output file')
    parser.add_argument('--score-cutoff',dest='score_cutoff', help='Cut-off for alignment score of leader (45)',default=45)
    parser.add_argument('--sample', help='sample id')
    parser.add_argument('--reads_file', help='the reads file from your 1st periscope run')

    args = parser.parse_args()
    print(args.score_cutoff)

    main(args)