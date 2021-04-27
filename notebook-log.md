# Botany 563 project notes

## Dataset
### Danionins
- The dataset for my project is mitochondrial DNA (mtDNA) from 10 species of Danionin fish (zebrafish and their relatives).
- These species include *Danio rerio*, *D*. *aesculapii*, *D*. *kyathit*, *D*. *tinwini*, *D*. *albolineatus*, *D*. *margaritatus*, *D*. *erythromicron*, *D*. *choprae*, *Microdevario kubotai*, and *Devario aequipinnatus*.
### Sequencing and annotation
- I sequenced *Danio rerio*, *D*. *aesculapii*, *D*. *kyathit*, *D*. *albolineatus*, *D*. *margaritatus*, and *Devario aequipinnatus* by using PCR to isolate mtDNA. UWBC performed library preparation, then used Illumina Hi-Seq 2500 or Nova-Seq to perform sequencing.
- Upon receipt of the paired-ends reads, I performed quality checks using FastQC, followed by adapter and quality trimming using fastp.
- I then assembled mitochondrial genomes out of the trimmed reads using SPAdes. I used SPAdes' read-depth filter to only get the reads that reached significant read depths (100,000 for Hi-Seq, 8,000 for Nova-Seq).
- I used MitoAnnotator to annotate genes from each genome then used a custom script to translate the amino acid sequences for each protein-coding gene.
### Other sequences
**From databases**
- Two species' mtDNA were downloaded from databases: *Danio erythromicron* (GenBank: AP011419.1) and *Microdevario kubotai* (RefSeq: NC_037360.1).
- MitoAnnotator was used to quickly annotate the genomes, and a custom script was used to translate each protein- coding gene.
#### From NCBI SRA projects
- Reads from whole-genome sequencing using Illumina NGS were downloaded from NCBI: *Danio tinwini* (SRA: ERX3311503, ERX3311504, ERX3311505, ERX3311506) and *D*. *choprae* (SRA: ERX3311487, ERX3311488, ERX3311489, ERX3311490).
- Reads were quality-checked with FastQC, and trimming for quality and adapters was performed using fastp.
- Reads files were concatenated, and iterative assembly for mitochondrial DNA was performed with MITObim using a zebrafish reference as bait.
- Assembled mtDNA was annotated using MitoAnnotator, and protein-coding genes were translated using a custom script.
### Plans for datasets
- Phylogenetic trees need to be created from each protein-coding gene, and possibly other conserved, or even non-coding regions.
- Each protein or codon sequence needs to be aligned, and proteins will be analyzed for selection. Sites of positive selection will be assessed as possible candidates for sites of functional difference between each gene in this phylogeny.



# Making a tree

### Extracted genes from MitoAnnotator files###
- H-strand genes (+), minus 12S rRNA, into plus.danio.fasta
`for i in AB aes kya nig tin alb mar ery cho kub son gnt dra sun mac esp`
`do`
`echo ">$i" >> plus.danio.fasta`
`cat mitoann.$i.genes.fa | sed -e "/^>/s/$/?/" -e "s/^>/@/" | tr -d "\n" | sed -e 's/?/\'$'\n/g' | sed -e 's/@/\'$'\n>/g' | sed -e '/>12S/,+1d' | sed -e '/>tRNA-Phe/,+1d' | grep -A 1 "+[\)]$" | sed -e '/--/,+0d' | grep -v "^>" | tr -d "\n" >> plus.danio.fasta`
`echo "" >> plus.danio.fasta`
`done`

### ClustalW2 alignment, fasta format
`~/software/clustalw2/clustalw2 -align -type=dna -outorder=input -output=fasta -infile=plus.danio.fasta -outfile=plus.danio.aln.fasta`

### ClustalW2 alignment, phy format
`~/software/clustalw2/clustalw2 -align -type=dna -outorder=input -output=phylip -infile=plus.danio.fasta -outfile=plus.danio.aln.phy`

### grep, mitoannotator annotations, and snapgene to create partitions
snapgene is essential in this part to match gaps through the find function.

* Genes by partition
    - tRNA (partition TRNA)
    trnaV - 1-72
    trnaL1 - 1806-1882
    trnaI - 2858-2929
    trnaM - 2930-2998
    trnaW - 4044-4117
    trnaD - 5671-5750
    trnaK - 6442-6517
    trnaG - 8152-8224
    trnaR - 8574-8644
    trnaH - 10324-10393
    trnaS2 - 10394-10463
    trnaL2 - 10464-10536
    trnaT - 13520-13593
    - 16S rRNA (partition RRNA)
    16S rRNA - 73-1805
    - Protein-coding genes by codon pos (Partitions POS1,POS2,POS3)
    nd1 - 1883-2857
    nd2 - 2999-4043
    cox1 - 4118-5670
    cox2 - 5751-6441
    atp8 - 6518-6682
    atp6 - 6683-7366
    cox3 - 7367-8151
    nd3 - 8225-8573
    nd4l - 8645-8941
    nd4 - 8942-10323
    nd5 - 10537-12375
    cytb - 12376-13519

- Example scripts (nd1 gene):
* Last nucleotide:
`cat plus.danio.aln.fasta | sed -e "/^>/s/$/?/" -e "s/^>/@/" | tr -d "\n" | sed -e 's/?/\'$'\n/g' | sed -e 's/@/\'$'\n>/g' | grep -A 1 "^>Danio_rerio" | grep -v "^>" | grep -o "[ACGT-]*ATATCGCCCTACCAATCGCACTAGCTGGTCTACCCCCACAAACATAA" | tr -d "\n" | wc -c`
2857
* First Nucleotide:
`cat plus.danio.aln.fasta | sed -e "/^>/s/$/?/" -e "s/^>/@/" | tr -d "\n" | sed -e 's/?/\'$'\n/g' | sed -e 's/@/\'$'\n>/g' | grep -A 1 "^>Danio_rerio" | grep -v "^>" | grep -o "[ACGT-]*ATGCTAGACATCTTAACAAGCCACTTAATTAACCCCCTAGCC" | sed s/"TGCTAGACATCTTAACAAGCCACTTAATTAACCCCCTAGCC"/""/g | tr -d "\n" | wc -c`
1883
- Final partitions:
DNA, RRNA=73-1805
DNA, TRNA=1-72,1806-1882,2858-2929,2930-2998,4044-4117,5671-5750,6442-6517,8152-8224,8574-8644,10324-10393,10394-10463,10464-10536,13520-13593
DNA, POS1=1883-2857/3,2999-4043/3,4118-5670/3,5751-6441/3,6518-6682/3,6683-7366/3,7367-8151/3,8225-8573/3,8645-8941/3,8942-10323/3,10537-12375/3,12376-13519/3
DNA, POS2=1884-2857/3,3000-4043/3,4119-5670/3,5752-6441/3,6519-6682/3,6684-7366/3,7368-8151/3,8226-8573/3,8646-8941/3,8943-10323/3,10538-12375/3,12377-13519/3
DNA, POS3=1885-2857/3,3001-4043/3,4120-5670/3,5753-6441/3,6520-6682/3,6685-7366/3,7369-8151/3,8227-8573/3,8647-8941/3,8944-10323/3,10539-12375/3,12378-13519/3

* Saved to file RT123.partition.txt

* Based off this example from Hal 1.3 (model will come later):
GTR+G+FO, NADH4=1-504/3,2-504/3
JC+I, tRNA=505-656
GTR+R4+FC, NADH5=657-898
HKY, NADH4p3=3-504/3

** Future trials may break further into resp complexes. Perhaps even single protein-coding genes.

### JModelTest2 on full alignment (no partitions)
- Looking for site model for plus.danio.aln.fasta
java -jar ~/software/jmodeltest-2.1.10/jModelTest.jar -d plus.danio.aln.fasta -s 11 -f -i -g 4 -AIC -BIC -DT -p -a -w > plus.danio.jmodel.txt

- Best model for AIC selected: 
* AIC     GTR+I+G         0.36    0.24    0.13    0.27    0.00    0.00      2.538  11.699   2.897   0.869  24.942   1.000    0.38    0.86
BIC     TIM2+I+G        0.36    0.24    0.13    0.27    0.00    0.00      2.914  12.528   2.914   1.000  26.610   1.000    0.38    0.86
DT      TIM2+I+G        0.36    0.24    0.13    0.27    0.00    0.00      2.914  12.528   2.914   1.000  26.610   1.000    0.38    0.86

### RAxML run for maximum-likelihood (no partitions)
`~/software/raxml-ng --msa plus.danio.aln.fasta --model GTR+I+G --prefix np25start --threads 32 --seed 2 --tree pars{25},rand{25}`

50 distinct starting trees!

### ModelTest-NG on partitioned data
`~/software/modeltest-ng-static −i plus.danio.aln.fasta -t ml -d nt -o plus.danio.partition.model.log -p 16 -q RT123.partition.txt`

------* The above did not work------------
modeltest-ng: Cannot parse the msa: 
              [900]: input file does not exist
modeltest-ng: You must specify an alignment file (-i)
Error: Invalid arguments
Try `modeltest-ng --help` for more information

### ModelTest-NG on partitioned data - absolute path
`~/software/modeltest-ng-static -i /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/plus.danio.aln.fasta -t ml -d nt -o /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/plus.danio.partition.model.log -p 16 -q /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/RT123.partition.txt`

Absolute paths do work. =)

* Partitions with models for RAXML. All based on best AIC.
GTR+I+G4, RRNA=73-1805
GTR+I+G4, TRNA=1-72,1806-1882,2858-2929,2930-2998,4044-4117,5671-5750,6442-6517,8152-8224,8574-8644,10324-10393,10394-10463,10464-10536,13520-13593
GTR+I+G4, POS1=1883-2857/3,2999-4043/3,4118-5670/3,5751-6441/3,6518-6682/3,6683-7366/3,7367-8151/3,8225-8573/3,8645-8941/3,8942-10323/3,10537-12375/3,12376-13519/3
TVM+I+G4, POS2=1884-2857/3,3000-4043/3,4119-5670/3,5752-6441/3,6519-6682/3,6684-7366/3,7368-8151/3,8226-8573/3,8646-8941/3,8943-10323/3,10538-12375/3,12377-13519/3
TIM2+I+G4, POS3=1885-2857/3,3001-4043/3,4120-5670/3,5753-6441/3,6520-6682/3,6685-7366/3,7369-8151/3,8227-8573/3,8647-8941/3,8944-10323/3,10539-12375/3,12378-13519/3

* Partitions saved into RT123.partition.model.txt

### RAxML using RT123.partition.model.txt partitioning scheme.
`~/software/raxml-ng --msa plus.danio.aln.fasta --model RT123.partition.model.txt --prefix part25start --threads 32 --seed 2 --tree pars{25},rand{25}`

### Trying new partitioning scheme with protein-coding genes divided by resp complex.

The below designations of ND, COX, or ATP do not refer to proteins, but define the respiratory complex the ND, COX, or ATP proteins are part of. The numbers define codon positions.

DNA, RRNA=73-1805
DNA, TRNA=1-72,1806-1882,2858-2929,2930-2998,4044-4117,5671-5750,6442-6517,8152-8224,8574-8644,10324-10393,10394-10463,10464-10536,13520-13593
DNA, ND1=1883-2857/3,2999-4043/3,8225-8573/3,8645-8941/3,8942-10323/3,10537-12375/3
DNA, ND2=1884-2857/3,3000-4043/3,8226-8573/3,8646-8941/3,8943-10323/3,10538-12375/3
DNA, ND3=1885-2857/3,3001-4043/3,8227-8573/3,8647-8941/3,8944-10323/3,10539-12375/3
DNA, CYB1=12376-13519/3
DNA, CYB2=12377-13519/3
DNA, CYB3=12378-13519/3
DNA, COX1=4118-5670/3,5751-6441/3,7367-8151/3
DNA, COX2=4119-5670/3,5752-6441/3,7368-8151/3
DNA, COX3=4120-5670/3,5753-6441/3,7369-8151/3
DNA, ATP1=6518-6682/3,6683-7366/3
DNA, ATP2=6519-6682/3,6684-7366/3
DNA, ATP3=6520-6682/3,6685-7366/3

* saved partitions into RTcomplex.partition.txt

### ModelTest-NG on partitioned + codon
`~/software/modeltest-ng-static -i /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/plus.danio.aln.fasta -t ml -d nt -o /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/plus.danio.complex.partition.model -p 16 -q /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/RTcomplex.partition.txt`

### Further partitioned EVERY GENE
DNA, tranV=1-72
DNA, trnaL1=1806-1882
DNA, trnaI=2858-2929
DNA, trnaM=2930-2998
DNA, trnaW=4044-4117
DNA, trnaD=5671-5750
DNA, trnaK=6442-6517
DNA, trnaG=8152-8224
DNA, trnaR=8574-8644
DNA, trnaH=10324-10393
DNA, trnaS2=10394-10463
DNA, trnaL2=10464-10536
DNA, trnaT=13520-13593
DNA, 16S=73-1805
DNA, nd1-1=1883-2857/3
DNA, nd1-2=1884-2857/3
DNA, nd1-3=1885-2857/3
DNA, nd2-1=2999-4043/3
DNA, nd2-2=3000-4043/3
DNA, nd2-3=3001-4043/3
DNA, cox1-1=4118-5670/3
DNA, cox1-2=4119-5670/3
DNA, cox1-3=4120-5670/3
DNA, cox2-1=5751-6441/3
DNA, cox2-2=5752-6441/3
DNA, cox2-3=5753-6441/3
DNA, atp8-1=6518-6682/3
DNA, atp8-2=6519-6682/3
DNA, atp8-3=6520-6682/3
DNA, atp6-1=6683-7366/3
DNA, atp6-2=6684-7366/3
DNA, atp6-3=6685-7366/3
DNA, cox3-1=7367-8151/3
DNA, cox3-2=7368-8151/3
DNA, cox3-3=7369-8151/3
DNA, nd3-1=8225-8573/3
DNA, nd3-2=8226-8573/3
DNA, nd3-3=8227-8573/3
DNA, nd4l-1=8645-8941/3
DNA, nd4l-2=8646-8941/3
DNA, nd4l-3=8647-8941/3
DNA, nd4-1=8942-10323/3
DNA, nd4-2=8943-10323/3
DNA, nd4-3=8944-10323/3
DNA, nd5-1=10537-12375/3
DNA, nd5-2=10538-12375/3
DNA, nd5-3=10539-12375/3
DNA, cytb-1=12376-13519/3
DNA, cytb-2=12377-13519/3
DNA, cytb-3=12378-13519/3

50 partitions!!!!!

* Saved partition scheme into gene.partition.txt

### ModelTest-NG on every gene
`~/software/modeltest-ng-static -i /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/plus.danio.aln.fasta -t ml -d nt -o /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/plus.danio.gene.partition.model -p 16 -q /mnt/sas0/AD/tjchamberlai/beast_stuff/aln/gene.partition.txt`

### Best AIC model for each in complex partition

GTR+G4, RRNA=73-1805
GTR+I+G4, TRNA=1-72,1806-1882,2858-2929,2930-2998,4044-4117,5671-5750,6442-6517,8152-8224,8574-8644,10324-10393,10394-10463,10464-10536,13520-13593
GTR+I+G4, ND1=1883-2857/3,2999-4043/3,8225-8573/3,8645-8941/3,8942-10323/3,10537-12375/3
TVM+I+G4, ND2=1884-2857/3,3000-4043/3,8226-8573/3,8646-8941/3,8943-10323/3,10538-12375/3
TIM1+I+G4, ND3=1885-2857/3,3001-4043/3,8227-8573/3,8647-8941/3,8944-10323/3,10539-12375/3
GTR+I+G4, CYB1=12376-13519/3
TIM1+I+G4, CYB2=12377-13519/3
TIM3+I+G4, CYB3=12378-13519/3
TPM3uf+G4, COX1=4118-5670/3,5751-6441/3,7367-8151/3
TPM3uf+I+G4, COX2=4119-5670/3,5752-6441/3,7368-8151/3
HKY+I+G4, COX3=4120-5670/3,5753-6441/3,7369-8151/3
TIM2ef+I+G4, ATP1=6518-6682/3,6683-7366/3
TrN+I, ATP2=6519-6682/3,6684-7366/3
TPM2uf+I+G4, ATP3=6520-6682/3,6685-7366/3

 * Saved into RTcomplex.partition.model.txt

### RAxML using RTcomplex.partition.model.txt partitioning scheme.
`~/software/raxml-ng --msa plus.danio.aln.fasta --model RTcomplex.partition.model.txt --prefix complex25 --threads 32 --seed 2 --tree pars{25},rand{25}`

### Best AIC model for each in EVERY GENE partition

Found each with grep:
`cat plus.danio.gene.partition.model.out | grep -A 2 "Best model according to AIC$" | sed -e '/Best/,+1d' | sed -e '/--/,+0d' | grep -o "[A-Za-z0-9+]*$"`

TPM2+G4, tranV=1-72
GTR+I+G4, trnaL1=1806-1882
TPM2+I, trnaI=2858-2929
GTR+G4, trnaM=2930-2998
TPM3uf+I+G4, trnaW=4044-4117
TPM2uf+I+G4, trnaD=5671-5750
TIM2ef+G4, trnaK=6442-6517
TPM2+G4, trnaG=8152-8224
TIM2+I+G4, trnaR=8574-8644
TIM3+I+G4, trnaH=10324-10393
TrN+I+G4, trnaS2=10394-10463
TPM3+G4, trnaL2=10464-10536
TIM3+G4, trnaT=13520-13593
HKY+I+G4, 16S=73-1805
TrN+I+G4, nd1-1=1883-2857/3
TIM2+G4, nd1-2=1884-2857/3
TIM3ef+I+G4, nd1-3=1885-2857/3
TPM2uf+I+G4, nd2-1=2999-4043/3
HKY+I+G4, nd2-2=3000-4043/3
TIM2+I, nd2-3=3001-4043/3
TPM3uf+I, cox1-1=4118-5670/3
TPM3uf+G4, cox1-2=4119-5670/3
HKY+I+G4, cox1-3=4120-5670/3
GTR+G4, cox2-1=5751-6441/3
HKY+I+G4, cox2-2=5752-6441/3
TrN+I+G4, cox2-3=5753-6441/3
TIM2ef+G4, atp8-1=6518-6682/3
TIM1+G4, atp8-2=6519-6682/3
TPM1uf+I+G4, atp8-3=6520-6682/3
TPM2uf+G4, atp6-1=6683-7366/3
SYM+I, atp6-2=6684-7366/3
TVM+G4, atp6-3=6685-7366/3
TIM1+I+G4, cox3-1=7367-8151/3
TPM2uf+G4, cox3-2=7368-8151/3
TIM2ef+G4, cox3-3=7369-8151/3
TPM1uf+G4, nd3-1=8225-8573/3
TrN+I+G4, nd3-2=8226-8573/3
GTR+I+G4, nd3-3=8227-8573/3
TIM2+I+G4, nd4l-1=8645-8941/3
TIM1+I+G4, nd4l-2=8646-8941/3
TPM1uf+G4, nd4l-3=8647-8941/3
HKY+G4, nd4-1=8942-10323/3
TPM1uf+I+G4, nd4-2=8943-10323/3
GTR+I+G4, nd4-3=8944-10323/3
TVM+G4, nd5-1=10537-12375/3
TrN+I+G4, nd5-2=10538-12375/3
TIM2ef+I+G4, nd5-3=10539-12375/3
TrN+G4, cytb-1=12376-13519/3
HKY+I+G4, cytb-2=12377-13519/3
TIM2ef+I+G4, cytb-3=12378-13519/3

* saved into gene.partition.model.txt

### RAxML using gene.partition.model.txt partitioning scheme.
`~/software/raxml-ng --msa plus.danio.aln.fasta --model gene.partition.model.txt --prefix gene25 --threads 32 --seed 2 --tree pars{25},rand{25}`

### RAxML bootstrap on no partitions.
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model GTR+I+G --prefix BSnp --seed 2 --threads 32`
* converged after 150 replicates

* manually increased replicates to 600
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model GTR+I+G --prefix BSnp600 --seed 2 --threads 16 --bs-trees 600`

* Assessing bootstrap convergence after inference - increased stringency.
`~/software/raxml-ng --bsconverge --bs-trees BSnp.raxml.bootstraps --prefix BSnp01 --seed 2 --threads 16 --bs-cutoff 0.01`
* no convergence after 150 trees.

* Assessing bootstrap convergence after inference - increased stringency - 600 trees.
`~/software/raxml-ng --bsconverge --bs-trees BSnp600.raxml.bootstraps --prefix BSnp01600 --seed 2 --threads 16 --bs-cutoff 0.01`

* No convergence after 600 trees.

* Performing bootstrap on an additional 600 trees. New seed value.
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model GTR+I+G --prefix BSnp600-2 --seed 333 --threads 16 --bs-trees 600`

* Concatenating trees. Assessing with increased stringency.
`cat BSnp600.raxml.bootstraps BSnp600-2.raxml.bootstraps > 1200BSnp.bootstraps`

`~/software/raxml-ng --bsconverge --bs-trees 1200BSnp.bootstraps --prefix BSnp01-1200 --seed 2 --threads 16 --bs-cutoff 0.01`
* Converged after 650 trees.

* Trying again with initial --bs-cutoff option
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model GTR+I+G --prefix BSnpinit01 --seed 2 --threads 32 --bs-cutoff 0.01`
* converged after 600 replicates

### RAxML bootstrap on RT123 partitions.
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model RT123.partition.model.txt --prefix BSRT123 --seed 2 --threads 32`
* Did not converge after 1000 replicates

### RAxML bootstrap on RTcomplex partitions.
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model RTcomplex.partition.model.txt --prefix BSRTcomplex --seed 2 --threads 16`
* Did not converge after 1000 replicates

### RAxML bootstrap on EVERY GENE partitions.
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model gene.partition.model.txt --prefix BSgene --seed 2 --threads 16`
* converged after 150 replicates

* Assessing after increased stringency of 0.01
`~/software/raxml-ng --bsconverge --bs-trees BSgene.raxml.bootstraps --prefix BSgene01 --seed 2 --threads 16 --bs-cutoff 0.01`

* WRF value decreasing rapidly. Will try with 350 more trees. New seed.
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model gene.partition.model.txt --prefix BSgene350 --seed 333 --threads 32 --bs-trees 350`

* Concatenated BSgene and BSgene350 trees and re-ran with 0.01 stringency.
`cat BSgene.raxml.bootstraps BSgene350.raxml.bootstraps > BSgene500.bootstraps`

`~/software/raxml-ng --bsconverge --bs-trees BSgene500.bootstraps --prefix BSgene500 --seed 2 --threads 32 --bs-cutoff 0.01`

* did not converge after 500 trees.

* 200 more trees
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model gene.partition.model.txt --prefix BSgene200 --seed 88 --threads 32 --bs-trees 200`

* Concatenating 200 trees to previous 500 trees, reassessing with 0.01 stringency
`cat BSgene200.raxml.bootstraps BSgene500.bootstraps > BSgene700.bootstraps`

`~/software/raxml-ng --bsconverge --bs-trees BSgene700.bootstraps --prefix BSgene700 --seed 2 --threads 32 --bs-cutoff 0.01`
* converged after 700 trees

* Experimenting with stringency cutoff in initial run. Tired of the hunt and peck.
`~/software/raxml-ng --bootstrap --msa plus.danio.aln.fasta --model gene.partition.model.txt --prefix BSgene01init --seed 2 --threads 16 --bs-cutoff 0.01`
* Converged after 800 replicates


* Trying the --all function with a --bs-cutoff 0.01 (new folder `all`)
`~/software/raxml-ng --all --msa plus.danio.aln.fasta --model gene.partition.model.txt --prefix gene-all --seed 2 --threads 32 --bs-metric fbp,tbe --bs-cutoff 0.01`
* REALLY EASY!!! All files in the folder were ready to go too!

### Mapping results onto trees
* starting with no partition, 2 trials
* in folder `np`
`~/software/raxml-ng --support --tree np25start.raxml.bestTree --bs-trees 1200BSnp.bootstraps --prefix 1200BSnpmap --threads 16`

`~/software/raxml-ng --support --tree np25start.raxml.bestTree --bs-trees BSnpinit01.raxml.bootstraps --prefix BSnpinit01map --threads 16`

* Now doing trees with EVERY GENE partition
* in folder `gene`
`~/software/raxml-ng --support --tree gene25.raxml.bestTree --bs-trees BSgene700.bootstraps --prefix BSgene700map --threads 16`

`~/software/raxml-ng --support --tree gene25.raxml.bestTree --bs-trees BSgene01init.raxml.bootstraps --prefix BSgene01initmap --threads 32`

## Trees and support look as good as can hope for using the EVERY GENE partition scheme!!!

# Running codeml to test for positive selection

### Some details about codeml and my running strategy
-codeml runs using a branch, or branches as forground, and keeps the other tree branches in the background. This way, it can look for selection in a single organism or clade. I'm looking across the Danionion lineage, one branch at a time, so I wrote a script named `every_branch.sh` that goes through the total number of branches (29 on my current tree), and assigns each of them as the forground one at a time. The foreground branch is defined by placing a #1 at the branch, which my script performs by looking for : and placing a space #1 at a branch before a run. On my laptop, I annotated the tree, gene25.raxml.bestTree, and use the image as a reference. The output for codeml is a long file in which I search for lnlog values and perform chidist tests using Excel. Significant p-values (I'm using <0.05) indicate statistically significant likelihoods, and I can report the amino acid positions on each branch that measured a high likelihood for positive selection (0.97 and above).

