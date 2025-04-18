# Whelpley_et_al_2018_SeaCucumberPhylo
Sea cucumber phylogenomics

## PLANNED ANALYSIS FOR TESTING HOLOTHURIAN PHYLOGENY AND DEVELOPING TARGET ENRICHMENT BAITS

Principle Investigators: Jessica Whelpley, Joseph Ryan, Gustav Paulay 
Draft or Version Number: v.2.2 
Date: 7 April 2025  
Note: this document will be updated (updates tracked through GitHub. See appendix)
 
## 1. INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_

Sea cucumbers are the most apomorphic echinoderms: bilaterally symmetrical worms with a reduced skeleton. They have evolved bizarre specializations such as anal suspension feeding, evisceration, sticky Cuvierian tubules that entangle attackers, and a “melting” body wall. They are abundant, ubiquitous in the benthos, intertidal to the deepest trenches from the poles to the equator, and include >1700 species in 25 families. They constitute the largest invertebrate fishery on Pacific islands, with stocks fully depleted throughout the tropics. A recent six gene molecular phylogeny showed substantial conflict with previous morphological-based relationships (Miller et al. 2017). Here, we present the pipeline of our phylogenomic analysis of Holothuroidea using 9 unpublished and 15 published holothurian transcriptomes and 16 published echinoderm transcriptomes that will serve as our outgroups.

### 1.2 _Rationale_ 

Inferring robust phylogenetic relationships is essential to understand how holothurians have adapted to the varied environments that they inhabit. Additionally, their evolutionary relationships can directly influence the development of fishing and conservation initiatives. A phylogeny built with hundreds of genes will help bring evolutionary resolution to this fascinating group of animals and provide an important framework to ask important questions related to their biodiversity and ecology. 

Part of this analysis will be to test the effects of data type, fossil choice, and models of molecular evolution on divergence dating. 

### 1.3 _Objectives_  

The overall objective is to test the higher-level relationships of holothurians recently proposed by Miller et al. 2017 by using more sequence data.

## 2. STUDY DESIGN AND ENDPOINTS 
#### 2.1 For this study we will use novel transcriptomes and publically available holothurian transcriptomes. 

2.1.1 We will use representatives from the other four classes of echinoderms (Asteroidea, Crinoidea, Echinoidea and Ophiuroidea) for our outgroups. Outgroup selection will be based on recent phylogenetic studies of Echinodermata as a whole group and studies of the individual classes. 

2.1.2 We will search for publically available transcriptomes from SRA database
(https://www.ncbi.nlm.nih.gov/sra) using the following parameters: 
```
taxa name [ORGN] AND Illumina BUTNOT genotyping BUTNOT metagenome BUTNOT miRNA
```

2.1.3 We will download and split the SRR sequences using SRA toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) then rename and compress the FASTQ files. (Note: `pigz` is a threaded version of `gzip`)

```
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR[number] > fqd.SRR[number].out 2> fqd.SRR[number].err &
```

```
mv SRR[number]_1.fastq SRR[number]_pass_1.fastq
```

```
pigz -9 SRR[number]_pass_1.fastq
```

#### 2.2 Transcriptome assembly and assessment. 

2.2.1 We will use BL-FILTER, part of the Agalma pipeline (Dunn et al. 2013), on the publically available transcriptomes and our transcriptomes to trim the adapters added during Illumina RNA-Seq. BL-FILTER and organizational steps for the downloaded transcriptomes are found below. Similar steps are applied to our transcriptomes, with slight variation at the beginning due to different file names.   

```
bl-filter-illumina -a -i ./SRR[number]_pass_1.fastq.gz -i ./SRR[number]_pass_2.fastq.gz -o SRR[number].1.fq -o SRR[number].2.fq -u SRR[number].unp.fq > blf.out 2> blf.err &
```

Copy the R1 reads to all_1.fq and then concatentate the unpaired reads to this file (this is how trinity allows for the incorporation of unpaired reads). Make a link from R2 file to all_2.fq so the file naming is consistent.

```
cp SRR[number].1.fq all_1.fq
```

```
cat SRR[number].unp.fq >> all_1.fq
```

```
ln -s SRR[number].2.fq.gz all_2.fq.gz
```

2.2.3 We will perform *de novo* transcriptome assembly with Trinity v2.4.0 (Grabherr et al. 2011).  

```
/usr/local/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 750G --CPU 10 --left ./SRR[number].1.fq.renamed --right ./SRR[number].2.fq.renamed --full_cleanup --normalize_reads --normalize_max_read_cov 30 > trin.out 2> trin.err &
```

#### 2.3 We will check assembly using BUSCO on the website gVolante (https://gvolante.riken.jp/) to assess transcriptome assembly. We will use the following parameters: 

Sequence type: coding/transcribed (nucleotide)
Choose an analysis program: BUSCO v2/v3
Ortholog set for BUSCO v2/v3: Eukaryota

#### 2.4 We will translate holothurian nucleotide transcriptome sequences into amino acid sequences with TransDecoder v3.0.1 (Haas & Papanicolaou 2012). We will set the –m flag to 50 and use the results from blast and hmmscan searches to inform the final TransDecoder prediction step.  

```
TransDecoder.LongOrfs -t [transcriptome_file] -m 50 > td.out 2> td.err
```

```
blastp -query longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > blastp.out 2> blastp.err 
```

```
hmmscan --cpu 1 --domtblout outfile.domtblout Pfam-A.hmm longest_orfs.pep > hs.out 2> hs.err
```

```
TransDecoder.Predict -t [transcriptome_file] --retain_pfam_hits out.domtblout --retain_blastp_hits out.blastp.out > tdp.out 2> tdp.err
```

#### 2.5 We will use the program alien_index (https://github.com/josephryan/alien_index) and version 0.02 of the alien_index database (a database of representative metazoan and non-metazoan sequences: http://ryanlab.whitney.ufl.edu/downloads/alien_index/) to identify and remove any potential contaminating/non-metazoan sequences. 

```
blastp -query [infile.pep.fa] -db ai.fa -outfmt 6 -max_target_seqs 1000 -seg yes -evalue 0.001 -out [file.out] > file.std 2> file.err
```

```
../alien_index --blast=[file_ai.out] --alien_pattern=ALIEN [out.alien_index] > ai.out 2> ai.err 
```

```
remove_aliens.pl [out.alien_index] [original_transcriptome.fa] > [filtered_transcriptome.fa] > ra.out 2> ra.err
```

#### 2.6 We will identify orthogroups across holothurian transcriptomes with OrthoFinder v2.1.2 (Emms & Kelly 2015).

```
orthofinder -f [dir_w_protein_fast_files] -op > of.out 2> of.err
```

##### Note 12/22/17: The first step generated 676 blastp files. We ran the following BLASTs and orthofinder command. 

```
blastp -outfmt 6 -evalue 0.001 -query [renamed_fasta_file_w_all_seqs] -db BlastDB -out outfile.txt > blastp.out 2> blastp.err
```

```
orthofinder -b [dir_w_blast_results] > ofb.out 2> ofb.err
```

#### 2.8 We will run the following commands to retain the orthogroup FASTA files that contain the user-specified minimum number of taxa (80%).

2.8.1 We will use Mafft v7.309 (Yamada et al. 2016) to align the sequences within each orthogroup. 

```mafft-linsi --localpair --maxiterate 1000 --thread 20 [infile] >mafft.out 2> mafft.err```

2.8.2 alignments will be refined using Gblockswrapper v0.03 (https://goo.gl/fDjan6)

```Gblockswrapper [infile.mafft] > outfile.mafft-gb > gbw.out 2> gbw.err```

2.8.3 Gblockswrapper sometimes leaves blank sequences that cause downstream issues. We removed any sequences with empty sequence lines. 

2.8.4 We will use IQ-TREE v1.6.3 (Nguyen et al. 2014) to estimate maximum-likelihood orthogroup gene trees. 

```iqtree -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix] > iq.out 2> iq.err```

2.8.5 We will use PhyloTreePruner v1.0 (Kocot et al. 2013) to prune the orthogroups with multiple sequences.

```java PhyloTreePruner [infile.tree] 25 [infile.align] 0.5 u > ptp.out 2> ptp.err```

2.8.6 We will use ```parapruner.pl``` to prune paraphyletic taxa from trees and alignments. This script can be found in the the scripts directory of this repository. In the tree directory we will run 

```ls -1 | perl -ne 'chomp; m/(OG[0-9]+)\.iq\.contree/; print "parapruner.pl --tree 01-TREES/$_ --aln 02-ALN/$1.mafft-gb.res --outdir $1 --delim=. > stdouterr/$1.out 2> stdouterr/$1.err &\n";' >> parapruner.sh```:q!

in order to generate the parapruner commands for each orthologous group. Prior to executing the parapruner shell script the deflines need to be changed from ```|``` to ```.```.

```perl -pi -e 's/\|/./g;' *.contree/.res```

2.8.7 We will use the genes that have been identified by ```parapruner.pl``` and re-run them through PhyloTreePruner in order to reduce the monophyletic duplicates to one copy. First we will copy the orthogroup directories from the paraprune directory and replace all of the ```.``` in the deflines of the parapruned mafft and contree files with ```|```.

```perl -pi -e 's/\./|/g;' ./OG*/*.mafft-gb.res.parapruned```

We will use ```post_parapruner.pl``` (scripts directory) to identify parapruned orthogroups that meet an 80% species cut-off (minumum 37 species) as well as ensuring that the sequence length >= 50, not including gaps.

```perl post_parapruner.pl 37 */*.mafft-gb.res.parapruned```

#### 2.9 We will concatenate single-copy loci filtered from step 2.7 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available in scripts directory of https://github.com/josephryan/JFR-PerlModules). We will edit the definition lines in each fasta file prior to running ```fasta2phylomatrix```.

```perl -pi.orig -e 's/\|.*$//;' *.fa```

```fasta2phylomatrix --dir=DIR_OF_FASTAFILES > matrix.fa```

#### 2.10 We will estimate species maximum-likelihood phylogeny.  

2.10.1 Concatenated matrix, Maximum Likelihood. We will estimate a bootstrapped (1000 ultrafast replicates) species phylogeny in IQ-TREE v1.6.3 using the concatenated dataset. We will use the flag -m TEST to find best partition scheme and the flag -nt AUTO to determine the appropriate number of threads to use for the analysis. The partition file will be created with the script ```fasta2phylomatrix```, which is available in the scripts repository.

```
iqtree –s [infile] –pre [outfile_prefix] –spp [partition file] -nt AUTO –m TEST –bb 1000 > iqo.out 2> iqo.err
```

## 3. WORK COMPLETED SO FAR WITH DATES  

13 September 2018 We have completed steps 2.1-2.9 prior to the release of phylotocol v 2.1

25 April 2018: We have completed steps 2.1-2.5 prior to the release of phylotocol v 2.0

5 October 2017: We have completed steps 2.1-2.2.1 prior to the release of phylotocol v 1.0 

22 December 2017: We have completed steps 2.1-2.6 prior to the release of phylotocol v 1.0


## 4. LITERATURE REFERENCED 

Dunn CW, Howison M, Zapata F. (2013). Agalma: an automated phylogenomics workflow. BMC Bioinformatics 14(1): 330. doi:10.1186/1471-2105-14-330

Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology, 16(1), 157.

Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PubMed PMID: 21572440.

Gblockswrapper: http://bit.ly/2svaKcR

Haas, B., & Papanicolaou, A. (2012). Transdecoder: https://transdecoder.github.io/

Kocot, K. M., Citarella, M. R., Moroz, L. L., & Halanych, K. M. (2013). PhyloTreePruner: a phylogenetic tree-based approach for selection of orthologous sequences for phylogenomics. Evolutionary Bioinformatics Online, 9, 429.

Miller, A. K., Kerr, A. M., Paulay, G., Reich, M., Wilson, N. G., Carvajal, J. I., & Rouse, G. W. (2017). Molecular phylogeny of extant Holothuroidea (Echinodermata). Molecular Phylogenetics and Evolution, 111, 110-131.

Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2014). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

Yamada, K. D., Tomii, K., & Katoh, K. (2016). Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. Bioinformatics, 32(21), 3246-3251.


## APPENDIX

Version&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;Date&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Significant Revisions

1.0 22 December 2017

2.0 25 April 2018

2.1 13 September 2018 Added sections 2.8.6 through 2.8.7 and associated perl scripts to the script repository

2.2 4 April 2025 – Removed several superfluous details. To simplifiy and publish this study, which was in danger of not ever being published, we have removed coalescant tree, ML with profile mixture model, dating, and ancestral state reconstruction analyses.

