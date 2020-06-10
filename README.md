# periscope

A tool to quantify sub-genomic RNA (sgRNA) expression in SARS-CoV-2 artic network amplicon sequencing data.
Initial classification of reads into sub-genomic or not based on https://www.biorxiv.org/content/10.1101/2020.04.10.029454v1.abstract

![alt text](https://github.com/sheffield-bioinformatics-core/periscope/blob/master/ocean.png "periscope")

# Requirements
periscope runs on MacOS and Linux. 


* conda
* Your raw fastq files from the artic protocol
* Install periscope

# Installation
```
git clone https://github.com/sheffield-bioinformatics-core/periscope.git && cd periscope
conda env create -f environment.yml
conda activate periscope
python setup.py install
```

# Execution
```
conda activate periscope

periscope \
    --fastq-dir <PATH_TO_DEMUXED_FASTQ> \
    --output-prefix <PATH_TO_OUTPUT> \
    --sample <SAMPLE_NAME> \
    --resources <PATH_TO_PERISCOPE_RESOURCES_FOLDER> \
    --score_cutoff <ALIGNMENT_CUTOFF_FOR_sgRNA> \
    --threads <THREADS_FOR_MAPPING>
```

# Pipeline overview
## Pre-Processing

* Collect demutiplexed pass fastqs
* Remap RAW artic protocol reads

## Counting

_This step takes roughly 1minute per 10k reads_
_Our median read count is ~250k and this will take around 25minutes_

* Read bam file
* Filter unmapped and secondary alignments
* Assign amplicon to read (using artic align_trim.py)
* Search for leader sequence
* Assign read to ORF
* Classify read (see Figure 1)
* Normalise a few ways

![alt text](https://github.com/sheffield-bioinformatics-core/periscope/blob/master/read_classification.png "periscope")<!-- .element height="10%" width="10%" -->
__Figure 1. Read Classification Algorithm__ 

## Normalisation

We have taken two approaches, a global normalisation based on mapped read counts or a local normalisation based on gRNA from the same amplicon.

* gRNA or sgRNA Per 100,000 mapped reads (gRPHT or sgRPHT)
    * We do this per amplicon and sum them in instances where multiple amplicons contribute to the final ORF count
* sgRNA can be normalised to per 1000 gRNA reads from the same amplicon - sgRPTg (normalising for amplicon efficiency differences)
    * There are things you need to note here:
        * multiple amplicons can contribute to reads which support the same sgRNA
        * we normalise on a per amplicon level and then sum these to get an overall normalised count


## Outputs:

#### a tab-delimited text file of 
- sample
- read id 
- read start position
- orf the read starts in
- pairwise alignment score
- read class (g/sgRNA) 
- amplicon number

#### a tab-delimited text file with each row an ORF an each column of
- raw count of putative sgRNA
- genomic reads for the respective amplicon
- total reads for the respective amplicon
- normalised count of putaive sgRNA (sgRNA/gRNA)

#### tagged bam file
- XS is the alignment score
- XA is the amplicon number
- XC is the assigned class (gDNA or sgDNA)
- XO is the orf assigned

Why periscope? SUB-genomic RNA, SUB-marine, periscope.
<div>Icons made by <a href="https://www.flaticon.com/authors/freepik" title="Freepik">Freepik</a> from <a href="https://www.flaticon.com/" title="Flaticon">www.flaticon.com</a></div>
