import argparse
import sys
import pysam
import os

def main(args):

    result={}

    counts = {
        "CTA": 0,
        "GAT": 0,
        "N501Y": {
            "A": 0,
            "T": 0
        },
        "coverage": {
            28279: 0,
            28280: 0,
            28281: 0,
            23062: 0,
        }
    }

    samfile = pysam.Samfile(args.bam, "rb")
    for pileupcolumn in samfile.pileup('MN908947.3', 28279, 28281):
        for pileupread in pileupcolumn.pileups:
            id = pileupread.alignment.query_name
            if pileupread.alignment.tags[12][1] != 'gRNA':
                continue
            if pileupread.query_position is None:
                continue
            if id not in result:
                # result[id]={}
                result[id]=[]
            if pileupcolumn.pos == 28279:
                # result[id][28279]=pileupread.alignment.query_sequence[pileupread.query_position]
                result[id].append(pileupread.alignment.query_sequence[pileupread.query_position])
            if pileupcolumn.pos == 28280:
                # result[id][28280] = pileupread.alignment.query_sequence[pileupread.query_position]
                result[id].append(pileupread.alignment.query_sequence[pileupread.query_position])
            if pileupcolumn.pos == 28281:
                # result[id][28281] = pileupread.alignment.query_sequence[pileupread.query_position]
                result[id].append(pileupread.alignment.query_sequence[pileupread.query_position])

    print(result)

    for pileupcolumn in samfile.pileup('MN908947.3', 23061, 23063):
        for pileupread in pileupcolumn.pileups:
            if pileupread.alignment.tags[12][1] != 'gRNA':
                continue
            if pileupread.query_position is None:
                continue
            if pileupcolumn.pos == 23062:
                # result[id][28279]=pileupread.alignment.query_sequence[pileupread.query_position]
                base=pileupread.alignment.query_sequence[pileupread.query_position]
                if base == "A":
                    counts["N501Y"]["A"]+=1
                if base == "T":
                    counts["N501Y"]["T"]+=1

    for read in result:
        if "".join(result[read]) == "CTA":
            counts["CTA"]+=1
        if "".join(result[read]) == "GAT":
            counts["GAT"]+=1

    sample=os.path.basename(args.bam).replace("_periscope.bam","")
    print("sample,N:D3:GAT,N:D3:CTA,S:N501:AAT,S:N501:TAT")
    print(sample+","+str(counts["GAT"])+","+str(counts["CTA"])+","+str(counts["N501Y"]["A"])+","+str(counts["N501Y"]["T"]))

    for i in result:
        if result[i] == ['C','T','A']:
            print(i)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='periscopre: Search for sgRNA reads in artic network SARS-CoV-2 sequencing data')
    parser.add_argument('--bam', help='bam file',default="The bam file of full artic reads")
    parser.add_argument('--threads', help='display progress bar', default=1)


    args = parser.parse_args()

    periscope = main(args)

    if periscope:
        print("all done", file=sys.stderr)
