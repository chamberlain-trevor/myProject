# Botany 563 project notes

## Dataset
### Danionins
- The dataset for my project is mitochondrial DNA (mtDNA)\ from 10 species of Danionin fish (zebrafish and their\ relatives).\
- These species include *Danio rerio*, *D*. *aesculapii*,\ *D*. *kyathit*, *D*. *tinwini*, *D*. *albolineatus*,\ *D*. *margaritatus*, *D*. *erythromicron*, *D*.\ *choprae*, *Microdevario kubotai*, and *Devario\ aequipinnatus*.\
### Sequencing and annotation
- I sequenced *Danio rerio*, *D*. *aesculapii*, *D*.\ *kyathit*, *D*. *albolineatus*, *D*. *margaritatus*, and\ *Devario aequipinnatus* by using PCR to isolate mtDNA.\ UWBC performed library preparation, then used Illumina\ Hi-Seq 2500 or Nova-Seq to perform sequencing.\
- Upon receipt of the paired-ends reads, I performed\ quality checks using FastQC, followed by adapter and\ quality trimming using fastp.\
- I then assembled mitochondrial genomes out of the\ trimmed reads using SPAdes. I used SPAdes' read-depth\ filter to only get the reads that reached significant\ read depths (100,000 for Hi-Seq, 8,000 for Nova-Seq).\
- I used MitoAnnotator to annotate genes from each\ genome\ then used a custom script to translate the amino\ acid sequences for each protein-coding gene.\
### Other sequences
**From databases**
- Two species' mtDNA were downloaded from databases:\ *Danio erythromicron* (GenBank: AP011419.1) and\ *Microdevario kubotai* (RefSeq: NC_037360.1).\
- MitoAnnotator was used to quickly annotate the genomes,\ and a custom script was used to translate each protein-\ coding gene.\
#### From NCBI SRA projects**
- Reads from whole-genome sequencing using Illumina NGS\ were downloaded from NCBI: *Danio tinwini* (SRA:\ ERX3311503, ERX3311504, ERX3311505, ERX3311506) and *D*.\ *choprae* (SRA: ERX3311487, ERX3311488, ERX3311489,\ ERX3311490).\
- Reads were quality-checked with FastQC, and trimming\ for quality and adapters was performed using fastp.\
- Reads files were concatenated, and iterative assembly\ for mitochondrial DNA was performed with MITObim using a\ zebrafish reference as bait.\
- Assembled mtDNA was annotated using MitoAnnotator, and\ protein-coding genes were translated using a custom script.\
### Plans for datasets
- Phylogenetic trees need to be created from each\ protein-coding gene, and possibly other conserved, or\ even non-coding regions.\
- Each protein or codon sequence needs to be aligned,\ and proteins will be analyzed for selection. Sites of\ positive selection will be assessed as possible candidates\ for sites of functional difference between\ each gene in this phylogeny.\
