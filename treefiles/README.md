### cuke_aln.fa.gz

Alignment used to generate tree (compressed)

### cuke_part.tre

Newick format of tree with bootstraps

### cuke_part_nex

Nexus formatted partition file

### IQ-TREE (v1.6.12) command used to generate tree

iqtree -s ./cuke_aln.fa -spp ./cuke_part_nex -nt AUTO -m TEST -bb 1000

