# Whelpley_et_al_2018_SeaCucumberPhylo
Sea cucumber phylogenomics

## PLANNED ANALYSIS FOR TESTING HOLOTHURIAN PHYLOGENY AND DEVELOPING TARGET ENRICHMENT BAITS

 Principle Investigator: Joseph Ryan, Gustav Paulay 
 Support Personnel: Jessica Whelpley  
 Draft or Version Number: v.1.0  
 Date: 22 December 2017  
 Note: this document will be updated (updates will be tracked through github)
 
## 1. INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_

Sea cucumbers are the most apomorphic echinoderms: bilaterally symmetrical worms with a reduced skeleton. They have evolved bizarre specializations such as anal suspension feeding, evisceration, sticky Cuvierian tubules that entangle attackers, and a “melting” body wall. They are abundant, ubiquitous in the benthos, from poles to equator, intertidal to the deepest trenches, and include >1700 species in 25 families. They are among the most conspicuous mobile invertebrates on reefs and the deep sea and constitute the largest invertebrate fishery on Pacific islands, with stocks fully depleted throughout the tropics. A recent six gene molecular phylogeny showed substantial conflict with previous morphological-based relationships (Miller et al. 2017). Here, we present a the pipeline of our phylogenomic analysis of Holothuroidea using 9 unpublished and 16 published transcriptomes.

### 1.2 _Rationale_ 

An in-depth understanding of phylogenetic relationships is essential to understand how holothurians have adapted to the varied environments that they inhabit. Additionally, comprehension of their evolutionary relationships can directly influence fishing and conservation initiatives. The combination of a backbone phylogeny built with hundreds of genes and high-quality baits for target enrichment, will help bring phylogenetic resolution to this fascinating group of animals and provide an important set of resources for systematists to conduct low-cost phylogenetic and population sampling. 

### 1.3 _Objectives_  

The overall objective is to test the recently proposed topology (Miller et al. 2017) of higher-level relationships of holothurians. With the sea cucumber orthogroups identified through our phylogenetic analysis, and our recently published draft *Australostichopus mollis* genome (Long et al. 2016) we will design a set of target-enrichment baits for Holothuroidea. 

## 2. STUDY DESIGN AND ENDPOINTS 

#### 2.1 Search for publically available transcriptomes from SRA database (https://www.ncbi.nlm.nih.gov/sra) using the following parameters: 
```
taxa name [ORGN] AND Illumina BUTNOT genotyping BUTNOT metagenome BUTNOT mRNA
```

2.1.1 Download and split the seventeen SRR sequences using SRA toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) then rename and compress the FASTQ files. 

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

2.2.1 Use BL-FILTER on the publically available transcriptomes and our transcriptomes, part of the Agalma pipeline (Dunn et al. 2013) to trim the adapters added during Illumina RNA-Seq. BL-FILTER and organizational steps for downloaded transcriptomes. Similar steps are applied to our transcriptomes, with slight variation at the beginning due to different file names.   

```
bl-filter-illumina -a -i ../../00-DATA/fastq/SRR[number]_pass_1.fastq.gz -i ../../00-DATA/fastq/[number]_pass_2.fastq.gz -o SRR[number].1.fq -o SRR[number].2.fq -u SRR[number].unp.fq > blf.out 2> blf.err &
```

Copy the R1 reads to all_1.fq and then concatentate the unpaired reads to this file (this is how trinity allows for the incorporation of unpaired reads). Make a link from R2 file to all_2.fq so the the file naming is consistent.

```
cp SRR[number].1.fq all_1.fq
```

```
cat SRR[number].unp.fq >> all_1.fq
```

```
ln -s SRR[number].2.fq.gz all_2.fq.gz
```

2.2.2 Use the script ```fix_names.pl``` (in the script repository) in order to fix the deflines of the SRA files so that they are formatted properly for transcriptome assembly with Trinity.

```
perl fix_names.pl SRR[number].1.fq.gz SRR[number].2.fq.gz SRR[number].unp.fq.gz > fix_sra_names.out 2> fix_sra_names.err
```

2.2.3 De novo transcriptome assembly with Trinity v2.4.0.  

```
/usr/local/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 750G --CPU 12 --left ../01-BL-FILTER/SRR[number].1.fq.renamed --right ../01-BL-FILTER/SRR[number].2.fq.renamed --full_cleanup --normalize_reads --normalize_max_read_cov 30 > trin.out 2> trin.err &
```
##### Note 12/22/17: Trinity runs were initially started with higher memory and central processors, however, we lowered those values due to server constraints.

2.2.4 Use RSEM (Li and Dewey, 2011) to measure the gene and isoform abundance. The scripts `align_and_estimate_abundance.pl` and `rsemgetbestseqs.py` are now included as a part of the Trinity package. 

```
align_and_estimate_abundance.pl --transcripts ../02-TRINITY/trinity_out_dir.Trinity.fasta --seqType fq --left [dir_with_BL-FILTER_results]/SRR[number].1.fq.gz --right [dir_with_BL-FILTER_results]/SRR[number].2.fq.gz --output_dir aea --est_method RSEM --aln_method bowtie2 --thread_count 100 --prep_reference > aea.out 2> aea.err &
```

```
rsemgetbestseqs.py ./aea/RSEM.isoforms.results [dir_with_Trinity_output]/trinity_out_dir.Trinity.fasta > rgbs.out 2> rgbs.err &
```

#### 2.3 Translate holothurian nucleotide transcriptome sequences into amino acid sequences with TransDecoder v3.0.0. We will set the –m flag to 50 and use the results from blast and hmmscan searches to inform the final TransDecoder prediction step.  

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

#### 2.4 The program [Alien Index](https://github.com/josephryan/alien_index) and a database of representative metazoan and non-metazoan sequences (http://ryanlab.whitney.ufl.edu/downloads/alien_index/) will allow us to remove any contaminating, non-metazoan sequences. 

```
blastp -query [infile.pep.fa] -db ai.fa -outfmt 6 -max_target_seqs 1000 -seg yes -evalue 0.001 -out [file.out] > file.std 2> file.err
```

```
../alien_index --blast=[file_ai.out] --alien_pattern=ALIEN [out.alien_index] > ai.out 2> ai.err 
```

```
perl /usr/local/bin/remove_aliens.pl [out.alien_index] [original_transcriptome.fa] > [filtered_transcriptome.fa] > ra.out 2> ra.err
```

#### 2.5 Identify orthogroups across holothurian transcriptomes with OrthoFinder v1.1.8. 
```
orthofinder -f [dir_w_protein_fast_files] -op > of.out 2> of.err
```
##### Note 12/22/17: The first step generated 676 blastp files. We used the script `blastp_parser.pl` [which can be found in the script repository] in order to separate the blastp files into smaller .sh files so they could be executed sequentially on the servers. 

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
##### and removed the duplicate transcriptome by placing a `#` in front of the duplicate species in the `SpeciesIDs.txt` file generated from the previous analysis. From this point forward, only 25 transcriptomes were utilized in the analysis.  

```
python trees_from_MSA.py [dir_w_orthofinder_results] > tfm.out 2> tfm.err
```

#### 2.6 Many of the transcriptomes had a large number of isoforms due to deep sequencing. Because of this, it is unlikely that we will have orthogroups with single-copy genes. We will use the file `Statistics_Overall_1.csv` generated by Orthofinder to determine a minimum-number of taxa that need to be present in each orthogroup. We will then examine the average number of sequences per species (ANSPS) in the orthogroups that contain the minimum-number of taxa present. We will then choose 1000 groups with ANSPS. 

##### Note 12/22/17: Based on the results of orthofinder, we determined the minimum-number of taxa to be 20 species, or 80% of the 25 species. 

#### 2.7 We will use the script ```filter_ogs_write_scripts.pl``` (available in the scripts directory of the DEEPC repository: https://github.com/josephryan/2017-DEEPC_Ctenophora/tree/master/scripts) to retain the orthogroup fasta files that contain the user-specified minimum number of taxa determined above. Lastly, ```filter_ogs_write_scripts.pl``` automates the following processes: 

2.7.1 sequences within each orthogroup are aligned using Mafft v7.309 

```mafft-linsi --localpair --maxiterate 1000 --thread 20 [infile] >mafft.out 2> mafft.err```

2.7.2 alignments are refined using Gblockswrapper v0.03 (https://goo.gl/fDjan6)

```Gblockswrapper [infile.mafft] > outfile.mafft-gb > gbw.out 2> gbw.err```

2.7.3 Gblockswrapper sometimes leaves blank sequences that cause downstream issues; the ```remove_empty_seqs``` script, available in the scripts directory of this repository, removes empty sequences and spaces from sequence lines. 

```remove_empty_seqs [outfile.mafft-gb] > res.out 2> res.err```

2.7.4 Maximum-likelihood orthogroup gene trees are estimated in IQTree v1.5.5 

```iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix] > iq.out 2> iq.err```

2.7.5 The orthogroups with multiple sequences will be pruned in PhyloTreePruner v1.0 

```java PhyloTreePruner [infile.tree] 25 [infile.align] 0.5 u > ptp.out 2> ptp.err```

##### Note 12/23/17: We plan on revisiting pruning methods by evaluating PhyloTreePruner and developing a custom pruning method after the preliminary results are presented at SICB 2018. 

#### 2.8 Concatenate (insert #) single-copy loci filtered from step 5 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available in the scripts directory of DEEPC: https://github.com/josephryan/2017-DEEPC_Ctenophora/tree/master/scripts). We will edit the definition lines in each fasta file (```perl -pi.orig -e 's/\|.*$//;' *.fa```) prior to running ```fasta2phylomatrix```.  

#### 2.9 Estimate species phylogeny using concatenated and coalescent gene tree/species tree methods.  

2.9.1 Concatenated matrix, Maximum Likelihood: estimate a bootstrapped (1000 ultrafast replicates) species phylogeny in IQtree v1.5.5 using the concatenated dataset. We will use the flag -m TEST to find best partition scheme and the flag -nt AUTO to determine the appropriate number of threads to use for the analysis. The partition file will be created with the script ```fasta2phylomatrix```, which is available in this respository.

```
iqtree-omp –s [infile] –pre [outfile_prefix] –spp [partition file] -nt AUTO –m TEST –bb 1000 > iqo.out 2> iqo.err
```

2.9.2 Coalescent-based phylogeny: estimate the species phylogeny using ASTRAL-II v4.11.1 and ASTRID v1.4. 

> i) Generate individual maximum-likelihood gene trees in IQtree. 

```
iqtree-omp –nt AUTO –s [infile] –pre [prefix_for_outfiles] –m MFP+MERGE –bb 1000 > iq.out 2> iq.err
```

> ii) ASTRAL-II constrains the search space to those species trees that derive their bipartitions from the input gene trees

```
java -Xmx1000000M -jar astral.jar -i [gene_trees_file] -o [output_file] > astral.out 2> astral.err
```

> iii) ASTRID uses a distance matrix generated from the input gene trees to estimate the species tree and is robust to missing data

```
ASTRID –i [infile] –o [outfile] –m bionj > astrid.out 2> astrid.err
```

> iv) Compute branch support using local posterior probabilities.  

#### 2.10 If there are conflicting species-tree topologies from 2., perform SOWH tests (implemented in sowhat v.0.36) to compare topologies. Any topologies that can be rejected with a P-Value <= 0.05 will be excluded from downstream analyes (but still reported in results). Constraint trees will be added to the phylotocol before running the tests.

2.10.1 example sowhat command line

```sowhat --constraint=[topology_to_be_tested] --aln=[alignment] --name=[name] --dir=[output_dir] --rax=[raxmlHPC-PTHREADS-SSE3 -T [num_threads]] ```


## 3. WORK COMPLETED SO FAR WITH DATES  

5 October 2017: We have completed steps 2.1-2.2.1 prior to the release of phylotocol v 1.0 

22 December 2017: We have completed steps 2.1-2.6 prior to the release of phylotocol v 1.0

## 4. LITERATURE REFERENCED 

Church, S. H., Ryan, J. F., & Dunn, C. W. (2015). Automation and evaluation of the SOWH Test with SOWHAT. Systematic Biology, 64(6), 1048-1058.

Dunn CW, Howison M, Zapata F. 2013. Agalma: an automated phylogenomics workflow. BMC Bioinformatics 14(1): 330. doi:10.1186/1471-2105-14-330

Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology, 16(1), 157.
Gblockswrapper: http://bit.ly/2svaKcR

Kocot, K. M., Citarella, M. R., Moroz, L. L., & Halanych, K. M. (2013). PhyloTreePruner: a phylogenetic tree-based approach for selection of orthologous sequences for phylogenomics. Evolutionary Bioinformatics Online, 9, 429.

Miller, A. K., Kerr, A. M., Paulay, G., Reich, M., Wilson, N. G., Carvajal, J. I., & Rouse, G. W. (2017). Molecular phylogeny of extant Holothuroidea (Echinodermata). Molecular Phylogenetics and Evolution, 111, 110-131.

Mirarab, S., & Warnow, T. (2015). ASTRAL-II: coalescent-based species tree estimation with many hundreds of taxa and thousands of genes. Bioinformatics, 31(12), i44-i52.

Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2014). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

Suyama, M., Torrents, D., & Bork, P. (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research, 34(suppl_2), W609-W612.

TransDecoder: https://transdecoder.github.io/

Vachaspati, P., & Warnow, T. (2015). ASTRID: accurate species trees from internode distances. BMC Genomics, 16(10), S3.

Yamada, K. D., Tomii, K., & Katoh, K. (2016). Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. Bioinformatics, 32(21), 3246-3251.

## APPENDIX

Version&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;Date&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Significant Revisions  
1.1  
1.2  
1.3  
1.4 
