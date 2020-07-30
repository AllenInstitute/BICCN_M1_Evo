#Steps to make TMM normalized BigWig files by Cluster

#Step 1: Make by cluster bam and CPM normalized bigWig
makeClusterBigWig_SNARE2.pl cluster_Barcodes.txt whitelist original_bam cluster_bw 2

#Step 2: make chr22_50bp_bins.bed
bedtools makewindows -w 50 -g refdata-cellranger-marmoset-3.2-GCF_000004665.1/Callithrix_jacchus.GCF_000004665.1/fasta/calJac3.fa.fai | grep ^chr22 > marmoset_chr22_50bp_bins.bed

#Step 3: Calculate EdgeR size factors
calcSizeFactors.r

#Step 4: Make TMM normalized bw files
sh makeBW.sh
