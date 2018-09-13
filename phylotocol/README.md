# Whelpley_et_al_2018_SeaCucumberPhylo
Sea cucumber phylogenomics

## PLANNED ANALYSIS FOR TESTING HOLOTHURIAN PHYLOGENY AND DEVELOPING TARGET ENRICHMENT BAITS

Principle Investigators: Jessica Whelpley, Joseph Ryan, Gustav Paulay 
Draft or Version Number: v.2.1 
Date: 13 September 2018  
Note: this document will be updated (updates will be tracked through GitHub)
 
## 1. INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_

Sea cucumbers are the most apomorphic echinoderms: bilaterally symmetrical worms with a reduced skeleton. They have evolved bizarre specializations such as anal suspension feeding, evisceration, sticky Cuvierian tubules that entangle attackers, and a “melting” body wall. They are abundant, ubiquitous in the benthos, intertidal to the deepest trenches from the poles to the equator, and include >1700 species in 25 families. They constitute the largest invertebrate fishery on Pacific islands, with stocks fully depleted throughout the tropics. A recent six gene molecular phylogeny showed substantial conflict with previous morphological-based relationships (Miller et al. 2017). Here, we present the pipeline of our phylogenomic analysis of Holothuroidea using 9 unpublished and 15 published holothurian transcriptomes and 16 published echinoderm transcriptomes that will serve as our outgroups.

### 1.2 _Rationale_ 

Inferring robust phylogenetic relationships is essential to understand how holothurians have adapted to the varied environments that they inhabit. Additionally, their evolutionary relationships can directly influence the development of fishing and conservation initiatives. A phylogeny built with hundreds of genes will help bring evolutionary resolution to this fascinating group of animals and provide an important framework to ask important questions related to their biodiversity and ecology. 

Part of this analysis will be to test the effects of data type, fossil choice, and models of molecular evolution on divergence dating. 

### 1.3 _Objectives_  

The overall objective is to test the higher-level relationships of holothurians recently proposed by Miller et al. 2017 by using more sequence data. Additionally, we will test hypotheses regarding ancestral state reconstruction and divergence dating.

## 2. STUDY DESIGN AND ENDPOINTS 
#### 2.1 For this study we will use novel transcriptomes and publically available holothurian transcriptomes. 

2.1.1 We will use representatives from the other four classes of echinoderms (Asteroidea, Crinoidea, Echinoidea and Ophiuroidea) for our outgroups. Outgroup selection will be based on recent phylogenetic studies of Echinodermata as a whole group and studies of the individual classes. Specific species transcriptome choice will be based on the 1.) phylogenetic position of the species (we would like to have a good representation of phylogenetic diversity), 2.) matrix occupancy and 3). fossil availability for the extant clade. 

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

2.2.2 Use the script ```fix_names.pl``` (available as a script within the JFR Perl Modules distribution (as of v1.1):  https://github.com/josephryan/JFR-PerlModules) to change the deflines of the SRA files so that they are formatted properly for transcriptome assembly with Trinity.

```
perl fix_names.pl SRR[number].1.fq.gz SRR[number].2.fq.gz SRR[number].unp.fq.gz > fix_sra_names.out 2> fix_sra_names.err
```

2.2.3 We will perform *de novo* transcriptome assembly with Trinity v2.4.0 (Grabherr et al. 2011).  

```
/usr/local/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 750G --CPU 10 --left ./SRR[number].1.fq.renamed --right ./SRR[number].2.fq.renamed --full_cleanup --normalize_reads --normalize_max_read_cov 30 > trin.out 2> trin.err &
```

2.2.4 Use RSEM (Li and Dewey, 2011) to measure the gene and isoform abundance. The script `align_and_estimate_abundance.pl` is included as a part of the Trinity package. The script `rsemgetbestseqs.py` is available from Warren Francis’ BitBucket page: https://bitbucket.org/wrf/sequences/src

```
align_and_estimate_abundance.pl --transcripts ./trinity_out_dir.Trinity.fasta --seqType fq --left ./SRR[number].1.fq.gz --right ./SRR[number].2.fq.gz --output_dir aea --est_method RSEM --aln_method bowtie2 --thread_count 100 --prep_reference > aea.out 2> aea.err &
```

```
rsemgetbestseqs.py ./aea/RSEM.isoforms.results ./trinity_out_dir.Trinity.fasta > rgbs.out 2> rgbs.err &
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

##### Note 12/22/17: The first step generated 676 blastp files. We used the script `blastp_parser.pl` (which can be found in the script repository for this directory) in order to separate the blastp files into smaller .sh files so they could be executed sequentially on the servers. 

```
blastp -outfmt 6 -evalue 0.001 -query [renamed_fasta_file_w_all_seqs] -db BlastDB -out outfile.txt > blastp.out 2> blastp.err
```

```
orthofinder -b [dir_w_blast_results] > ofb.out 2> ofb.err
```

##### Note: 12/22/17: At this point we recognized an incorrectly labeled transcriptome from the NCBI database indicating that there was a duplicate transcriptome. I re-ran: 

```
orthofinder -b [dir_w_blast_results] > ofb.out 2> ofb.err
```

##### and removed the duplicate transcriptome by placing a `#` in front of the duplicate species in the `SpeciesIDs.txt` file generated from the previous analysis. From this point forward, only 24 transcriptomes were used in the analysis.  

```
python trees_from_MSA.py [dir_w_orthofinder_results] > tfm.out 2> tfm.err
```

#### 2.7 Due to the depth of sequencing, many of the transcriptomes will have a large number of isoforms. Because of this, it is unlikely that we will have orthogroups with single-copy genes. We will use the file `Statistics_Overall_1.csv` generated by Orthofinder to determine a minimum-number of taxa that need to be present in each orthogroup. The minimum-number of taxa we plan to use will include ~80% of all taxa in the study.

#### 2.8 We will use the script ```filter_ogs_write_scripts.pl``` (available in the scripts directory of this repository) to retain the orthogroup FASTA files that contain the user-specified minimum number of taxa determined above. Lastly, 

```filter_ogs_write_scripts.pl``` automates the following processes: 

2.8.1 We will use Mafft v7.309 (Yamada et al. 2016) to align the sequences within each orthogroup. 

```mafft-linsi --localpair --maxiterate 1000 --thread 20 [infile] >mafft.out 2> mafft.err```

2.8.2 alignments will be refined using Gblockswrapper v0.03 (https://goo.gl/fDjan6)

```Gblockswrapper [infile.mafft] > outfile.mafft-gb > gbw.out 2> gbw.err```

2.8.3 Gblockswrapper sometimes leaves blank sequences that cause downstream issues; the ```remove_empty_seqs``` script (available in the scripts directory of this repository) removes empty sequences and spaces from sequence lines. 

```remove_empty_seqs [outfile.mafft-gb] > res.out 2> res.err```

2.8.4 We will use IQ-TREE v1.6.3 (Nguyen et al. 2014) to estimate maximum-likelihood orthogroup gene trees. 

```IQ-TREE-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix] > iq.out 2> iq.err```

2.8.5 We will use PhyloTreePruner v1.0 (Kocot et al. 2013) to prune the orthogroups with multiple sequences.

```java PhyloTreePruner [infile.tree] 25 [infile.align] 0.5 u > ptp.out 2> ptp.err```

2.8.6 After the commands generated from ```filter_ogs_write_scripts.pl``` have been executed, we will copy the output IQ-tree ```.contree``` files and the alignment files ```.mafft-gb.res``` into a new, separate directories and then run script ```parapruner.pl```. This script ignores monophyletic duplicates in the gene tree and identifies and if there are more than a set number of paraphyletic duplicates, the script will remove them. This script can be found in the the scripts repository of this directory. In the tree directory we will run 

```ls -1 | perl -ne 'chomp; m/(OG[0-9]+)\.iq\.contree/; print "perl /usr/local/bin/parapruner.pl --tree 01-TREES/$_ --aln 02-ALN/$1.mafft-gb.res --outdir $1 --delim=. > stdouterr/$1.out 2> stdouterr/$1.err &\n";' >> parapruner.sh```:q!

in order to generate the parapruner commands for each orthologous group. Prior to executing the parapruner shell script the deflines need to be changed from ```|``` to ```.```.

```perl -pi -e 's/\|/./g;' *.contree/.res```

2.8.7 We will use the genes that have been identified by ```parapruner.pl``` and re-run them through PhyloTreePruner in order to reduce the monophyletic duplicates to one copy. First we will copy the orthogroup directories from the paraprune directory and replace all of the ```.``` in the deflines of the parapruned mafft and contree files with ```|```.

```perl -pi -e 's/\./|/g;' ./OG*/*.mafft-gb.res.parapruned```

Then we will execute the ```post_parapruner.pl``` script which can be found in the scripts repository of this directory. This script will identify the parapruned orthogroups that meet the 80% species cut-off value we previously used in 2.7 (set when writing the command) as well as ensuring that the sequence length meets a specified length (set by modifying the perl script), not including gaps. For this project we will set our sequence length to 50. 

```perl post_parapruner.pl [species cut-off value] */*.mafft-gb.res.parapruned```

#### 2.9 We will concatenate the single-copy loci filtered from step 2.7 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available as a script within the JFR Perl Modules distribution (as of v1.1):  https://github.com/josephryan/JFR-PerlModules) . We will edit the definition lines in each fasta file (```perl -pi.orig -e 's/\|.*$//;' *.fa```) prior to running ```fasta2phylomatrix```. 

#### 2.10 We will estimate species phylogeny using concatenated maximum-likelihood (ML) and coalescent gene tree/species tree methods. We will consider our maximum-likelihood tree as our working hypothesis, and we will use our coalescent consensus tree to assess the robustness of our ML topology. The differences and similarities between the trees will be discussed.  

2.10.1 Concatenated matrix, Maximum Likelihood. We will estimate a bootstrapped (1000 ultrafast replicates) species phylogeny in IQ-TREE v1.6.3 using the concatenated dataset. We will use the flag -m TEST to find best partition scheme and the flag -nt AUTO to determine the appropriate number of threads to use for the analysis. The partition file will be created with the script ```fasta2phylomatrix```, which is available in the scripts repository.

```
IQ-TREE-omp –s [infile] –pre [outfile_prefix] –spp [partition file] -nt AUTO –m TEST –bb 1000 > iqo.out 2> iqo.err
```
2.10.2 As an alternative method, we will apply the posterior mean site frequency (PMSF) model in IQ-TREE (Wang et al. 2017) to this dataset. We will use the ML tree generated from 2.10.1 as our guide tree for this analysis. 

```
iqtree -s <alignment> -m LG+C60+F+G -ft <guide_tree> -b100
```

2.10.3 Coalescent-based phylogeny. We will estimate the species phylogeny using ASTRAL-III (Zhang et al. 2017) v5.1.1 and ASTRID v1.4 (Vachaspati & Warnow 2015). 

> i) Generate individual maximum-likelihood gene trees in IQ-TREE. 

```
IQ-TREE-omp –nt AUTO –s [infile] –pre [prefix_for_outfiles] –m MFP+MERGE –bb 1000 > iq.out 2> iq.err
```

> ii) ASTRAL-III constrains the search space to those species trees that derive their bipartitions from the input gene trees

```
java -Xmx1000000M -jar astral.jar -i [gene_trees_file] -o [output_file] > astral.out 2> astral.err
```

> iii) ASTRID uses a distance matrix generated from the input gene trees to estimate the species tree and is robust to missing data.

```
ASTRID –i [infile] –o [outfile] –m bionj > astrid.out 2> astrid.err
```

> iv) Compute branch support using local posterior probabilities.  

#### 2.11 If there are conflicting species-tree topologies from 2.10, we will perform AU tests (implemented in sowhat v0.36) to compare topologies. Any topologies that can be rejected with a P-Value <= 0.05 will be excluded from downstream analyses (but still reported in results). Constraint trees will be added to the phylotocol before running the tests.

#### 2.12 Divergence dating. In this study we will explore the effects of tree topology (based on different data sources), fossil choice, and model choice on divergence estimates. Divergence will be estimated with BEAST v2.4.8 (Bouckaert et al. 2014)  using our maximum-likelihood tree (from 2.10.1). BEAUti will be used to generate the necessary .xml file for analysis. We will conduct the same divergence dating analyses performed in Miller et al. (2017) with some variation. Miller et al. (2017) constrained specific nodes and allowed BEAST to infer the remaining relationships during the analysis. Unlike their paper, we will use our ML tree generated from greatly increased sequence sampling and we will not allow BEAST to infer phylogenetic relationships. In addition, we will vary two other parameters, fossil inclusion and model choice, to see their effects on these analyses. 

2.12.1 The following describes the four BEAST analyses that will be run. All analyses will be run on our ML topology (from 2.10.1):

2.12.1.A: include the same set of fossils and parameters from Miller et al. (2017).

2.12.1.B: include the same set of fossils from the Miller et al. (2017) with parameters from our IQ-Tree analysis (see below). 

2.12.1.C: include a curated subset of fossils from Miller et al. (2017) with the original parameters from the Miller study. 

2.12.1.D: include a curated subset of fossils from Miller et al. (2017) with parameters from our IQ-Tree analysis (see below). 

Miller et al. (2017) used an uncorrelated log normal rate change, GTR + I + Γ model of substitution (partitioned) and Yule model of divergences. We will use similar parameters (uncorrelated relaxed clock with a lognormal distribution, Yule model of divergence) but will use our partitioned models of evolution from our IQ-TREE analysis (section 2.10.1). Initial analyses will be run for a minimum of 10 million generations and chain sampling will occur every 1,000 generations. We will consider ESS values of 200 or greater as an indication of convergence. We will initially run 10 million generations; if our ESS values are below 200, we will continue the analysis until it reaches 50 million generations. If our ESS values are still above 200; we will continue the analysis until it reaches 100 million generations. If after 100 million generations, ESS values have not reached 200, we will determine that the analyses have not converged and report this as the main finding regarding these analyses. The ESS values will be visualized in Tracer v1.4 (Rambaut & Drummond 2007) to observe information regarding the posterior, prior, likelihood and the parameters that we are estimating. 

2.12.2 Statistical comparisons of BEAST analyses. 

HYPOTHESIS 1: The comparison of p-values (derived from the comparison of nodes) will be significant for all three tests. The distribution of ages at shared nodes between compared analyses (see below) will be significantly different from each other at most nodes.

HYPOTHESIS 2: Variations in fossils included for fossil calibration (comparing 2.12.1.A and 2.12.1.C) will have the most significant effect on date discrepancies (when comparing average p-values between each comparison below).

To test the effects of topology, fossil choice, and model on molecular dating, we will compare date ranges of conserved ancestral nodes between the four BEAST analyses (A-D) proposed above. To test hypotheses 1 and 2, we will conduct two-tailed unpaired T-tests on each conserved ancestral node between all of our BEAST analyses (i.e. A:B, A:C, A:D, B:C, B:D, C:D). If there is a significant difference between the dates at any of the nodes, drawing a conclusion from these analyses for that node is problematic. If the results between the various analyses are the same at each node, the analysis and the data are robust and a conclusion can be made. These results will be discussed in detail as part of the supplementary information for the paper. 

Additionally, we will conduct two-tailed unpaired T-tests between the nodes that are shared in the Miller et al. (2017) BEAST analysis and our new BEAST analyses to rank the influence of the parameters we have decided to use for this study. To test the topology robustness we will compare the conserved nodes in the Miller et al. (2017) tree to our BEAST analysis A. We will test model robustness by comparing the conserved nodes from analysis A to B and analysis C to D. For fossil robustness, we will compare the conserved nodes from analysis A to C and B to D. We will set the topology, model of evolution and fossil as a variable (e.g. x = topology robustness, y = model of evolution robustness, z = fossil choice robustness) and use the T-test values to determine the rank of these variables. For example, if z is the highest ranked variable, then fossil choice has the greatest influence on the result of the analysis. If topology (x) is ranked the highest or the lowest, we need to perform an additional ranking of the fossil (z) and model of evolution (y). 

#### 2.13 Ancestral state reconstruction. We will use methodology similar to that found in Sasson et al. 2017. For completed R script using their methodology, visit their GitHub: https://github.com/josephryan/2017b_Sasson_and_Ryan. 

We will conduct ancestral state reconstructions on both the phylogram generated in section 2.10.1 and the ultrametric tree generated from our BEAST analysis. The differences between the results from these two trees will be discussed in the paper. We will use stochastic character mapping in Phytools v0.6-44 (Revell 2012) to estimate ancestral character states of holothurians using discrete variables. In order to run Phytools, we need to determine the model necessary to perform a likelihood ratio test. In R v3.5.0 we will use the ace function in ape v5.0 to estimate ancestral character states using a one-parameter equal rates model (ER), a three-parameter symmetrical model (SYM), and a six-parameter all-rates-different model (ARD) and perform a likelihood ratio test on the results. 

Using stochastic character mapping in Phytools, we will use the function make.simmap in order to calculate the likelihood for the 11 discrete variables we are testing. 

The likelihood ratio test and simmmap scripts can be found in the script repository for this project. 

#### 2.14 Proposed figures and tables for paper.
 
Figure 1: Images of some of the sea cucumbers used for this analysis.
Figure 2: Concatenated maximum-likelihood tree with support values from each 3 analyses (ML, C60, coalescent) 
Figure 3: Chart comparing BEAST analyses (BEAST trees for each analysis will be included in the supplementary material). 
Figure 4: Ancestral-state reconstruction of 11 discrete variables (phylogram and ultrametric trees).  
Table 1: Voucher information and accession numbers for specimens.
Table 2: Results of AU test (if performed). 
Table 3a: Result of T-tests comparing conserved nodes from BEAST analysis.
Table 3b: Result of T-tests comparing conserved nodes from the Miller et al. (2017) analysis and our analyses. 
Supplemental Figure 1: IQ-TREE C60 tree 
Supplemental Figure 2: Coalescent based tree.

## 3. WORK COMPLETED SO FAR WITH DATES  

13 September 2018 We have completed steps 2.1-2.9 prior to the release of phylotocol v 2.1

25 April 2018: We have completed steps 2.1-2.5 prior to the release of phylotocol v 2.0

5 October 2017: We have completed steps 2.1-2.2.1 prior to the release of phylotocol v 1.0 

22 December 2017: We have completed steps 2.1-2.6 prior to the release of phylotocol v 1.0

UPDATE PRIOR TO RE-PUBLISHING

## 4. LITERATURE REFERENCED 
Bouckaert R, Heled J, Kühnert D, Vaughan T, Wu C-H, Xie D, et al. (2014) BEAST 2: A Software Platform for Bayesian Evolutionary Analysis. PLoS Comput Biol 10(4): e1003537.

Church, S. H., Ryan, J. F., & Dunn, C. W. (2015). Automation and evaluation of the SOWH Test with SOWHAT. Systematic Biology, 64(6), 1048-1058.

Dunn CW, Howison M, Zapata F. (2013). Agalma: an automated phylogenomics workflow. BMC Bioinformatics 14(1): 330. doi:10.1186/1471-2105-14-330

Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology, 16(1), 157.

Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PubMed PMID: 21572440.

Gblockswrapper: http://bit.ly/2svaKcR

Huelsenbeck JP, Nielsen R, Bollback JP. Stochastic mapping of morphological characters. Sys Biol. 2003;52(2):131–58.
59.

Kocot, K. M., Citarella, M. R., Moroz, L. L., & Halanych, K. M. (2013). PhyloTreePruner: a phylogenetic tree-based approach for selection of orthologous sequences for phylogenomics. Evolutionary Bioinformatics Online, 9, 429.

Miller, A. K., Kerr, A. M., Paulay, G., Reich, M., Wilson, N. G., Carvajal, J. I., & Rouse, G. W. (2017). Molecular phylogeny of extant Holothuroidea (Echinodermata). Molecular Phylogenetics and Evolution, 111, 110-131.

Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2014). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

Paradis E, Claude J, Strimmer K. APE: analyses of phylogenetics and evolution in R language. Bioinformatics. 2004;20(2):289–90.

Rambaut A, & Drummond AJ (2007) Tracer v1.4, Available from http://beast.bio.ed.ac.uk/Tracer

Revell LJ. phytools: an R package for phylogenetic comparative biology (and other things). Methods Ecol Evol. 2012;3(2):217–23.
61.

Sasson, D. A., & Ryan, J. F. (2017). A reconstruction of sexual modes throughout animal evolution. BMC evolutionary biology, 17(1), 242

Suyama, M., Torrents, D., & Bork, P. (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research, 34(suppl_2), W609-W612.

Haas, B., & Papanicolaou, A. (2012). Transdecoder: https://transdecoder.github.io/

Vachaspati, P., & Warnow, T. (2015). ASTRID: accurate species trees from internode distances. BMC Genomics, 16(10), S3.

Yamada, K. D., Tomii, K., & Katoh, K. (2016). Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. Bioinformatics, 32(21), 3246-3251.

Zhang, Chao, Erfan Sayyari, and Siavash Mirarab. “ASTRAL-III: Increased Scalability and Impacts of Contracting Low Support Branches.” In Comparative Genomics: 15th International Workshop, RECOMB CG, 2017. 


## APPENDIX

Version&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;Date&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Significant Revisions  
1.0 22 December 2017
2.0 25 April 2018
2.1 13 September 2018 


