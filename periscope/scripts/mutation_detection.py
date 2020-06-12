import vcf
import pysam
import argparse
import pandas as pd
from plotnine import *




def check_position_in_bam(bam,position):

    result={}
    for pileupcolumn in bam.pileup("MN908947.3",position,position+1):
        # print("\ncoverage at base %s = %s" %
        #       (pileupcolumn.pos, pileupcolumn.n))
        if pileupcolumn.pos==position-1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    read_class = pileupread.alignment.get_tag("XC")
                    # print(read_class)
                    if read_class not in result:
                        result[read_class] = {}

                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base not in result[read_class]:
                        result[read_class][base]=0
                    result[read_class][base] += 1


    return result

def main(args):
    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    bam = pysam.AlignmentFile(args.bam, "rb")
    result = []
    poses = []
    for record in vcf_reader:
        check = check_position_in_bam(bam, record.POS)
        poses.append(record.POS)
        result.append(pd.DataFrame.from_dict(check))
    df = pd.concat(result,keys=poses)
    print(df)

    df.reset_index(inplace = True)
    df['level_0']=df['level_0'].astype(str)
    print(df)

    df = pd.melt(df,id_vars=['level_0','level_1'])

    pos_list = df['level_0'].value_counts().index.tolist()
    pos_cat = pd.Categorical(df['level_0'], categories=pos_list)
    df = df.assign(pos_cat=pos_cat)
    x=(ggplot(df, aes(x='pos_cat',y='value',color='level_1',fill='level_1'))) + \
      geom_bar(stat='identity',position='fill') + \
      facet_wrap(['variable'],ncol=1)

    x.save(filename=args.output_prefix+"_base_counts.png", height=5, width=5, units='in', dpi=300)

    df = df.rename(columns={"level_0": "position", "level_1": "base"})
    df = df.assign(sample=args.sample)
    print(df)
    df.to_csv(args.output_prefix+"_base_counts.csv",index=False)

    print(x)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='periscopre: Search for sgRNA reads in artic network SARS-CoV-2 sequencing data')
    parser.add_argument('--periscope-bam',dest="bam", help='bam file',default="/tmp/SHEF-C0C35_periscope.bam")
    parser.add_argument('--output-prefix',dest='output_prefix', help='Prefix of the output file',default="/tmp/TEST")
    parser.add_argument('--vcf', dest='vcf', help='Cut-off for alignment score of leader (45)',default="/tmp/SHEF-C0C35.pass.vcf")
    parser.add_argument('--sample', help='sample id',default="TEST")

    args = parser.parse_args()
    main(args)