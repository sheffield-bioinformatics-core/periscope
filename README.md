# periscope

A tool to quantify sub-genomic RNA expression in SARS-CoV-2 artic network amplicon sequencing data.

![alt text](https://github.com/sheffield-bioinformatics-core/periscope/blob/master/ocean.png "periscope")


Pipeline overview
* Read bam file
* Assign amplicon to read (using artic align_trim.py)
* Search for leader sequence in read to classify into sgRNA or gRNA
* Classify sgRNAs into ORFs based on start position
* Classy which amplicon an ORF can come from
* Normalise based on total read count per amplicon



Why periscopre? SUB-genomic RNA, SUB-marine, periscope.
<div>Icons made by <a href="https://www.flaticon.com/authors/freepik" title="Freepik">Freepik</a> from <a href="https://www.flaticon.com/" title="Flaticon">www.flaticon.com</a></div>
