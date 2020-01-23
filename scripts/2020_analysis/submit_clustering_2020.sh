#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=2G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace
# #$ -l h=blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# MergedSeqs=$OutDir/${Locus}_concatenated.temp.fa
# OutDir=OutDir=clustering/$Locus
# Prefix=$Locus
MergedSeqs=$1
OutDir=$2
Prefix=$3

CurDir=$PWD
WorkDir=${TMPDIR}/cluster_reads
mkdir -p $WorkDir
cd $WorkDir



# The sequences are reduced to only unique sequences ("dereplicated") to improve computation time:
vsearch --derep_fulllength $CurDir/$MergedSeqs --output ${Prefix}_derep.fa --sizeout --minseqlength 50
# Chimeras are detected in a de novo fashion, and excluded:
# vsearch --uchime_denovo ${Prefix}_derep.fa --nonchimeras ${Prefix}_nochimeras.fa
# Sequences are sorted by size, and singletons are removed (they are mapped back on later):
# vsearch --sortbysize ${Prefix}_nochimeras.fa --output ${Prefix}_sorted.fa --minsize 2
vsearch --sortbysize ${Prefix}_derep.fa --output ${Prefix}_sorted.fa --minsize 2
# Operational taxonomic units are clustered at 97% sequence identity, and the consensus sequences for the OTUs are output to a fa file:
vsearch --cluster_smallmem ${Prefix}_sorted.fa --id 0.97 --consout ${Prefix}_cluster_vsearch_identity_unparsed.fa --usersort --relabel OTU
# The OTU names contain a wealth of information that is not required for downstream processing. To reduce complications, we can simplify the names with awk:
awk 'BEGIN{OFS="";ORS="";count=0}{if ($0~/>/){if (NR>1) {print "\n"} print ">" count "\n"; count+=1;} else {print $0;}}' ${Prefix}_cluster_vsearch_identity_unparsed.fa > ${Prefix}_cluster_vsearch_OTUs.fa

vsearch --cluster_unoise ${Prefix}_sorted.fa --minsize 8 --unoise_alpha 2.0 --consout ${Prefix}_cluster_vsearch_zOTUs_unparsed.fa

vsearch --uchime3_denovo ${Prefix}_cluster_vsearch_zOTUs_unparsed.fa --nonchimeras ${Prefix}_cluster_vsearch_zOTUs.fa --relabel OTU


# vsearch --usearch_global stool_sequences.fa --db cluster/rep_set_relabel.fa --strand both --id 0.97 --uc cluster/map.uc --threads 2
# wget https://github.com/neufeld/MESaS/raw/master/scripts/mesas-uc2clust
# python mesas-uc2clust cluster/map.uc cluster/seq_otus.txt

#### Clustering (Cluster dereplicated seqeunces and produce OTU fa (also filters for chimeras))

# ---
# Identify OTUs
# ---
# With 97% clustering, an OTU sequence should be at least 3% different from all
# other OTUs, and should be the most abundant sequences in its neighborhood.
# This is done by the cluster_otus command, which is an implementation of the
# UPARSE algorithm.
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -cluster_otus ${Prefix}_derep.fa -otus ${Prefix}_cluster_usearch_OTUs.fa -relabel OTU -minsize 4


# ---
# Identify zOTUs
# ---
# Denoising attempts to identify all correct biological sequences in the reads.
# This is done by the unoise3 command, which is an implementation of the UNOISE
# algorithm. A denoised sequence is called a "ZOTU" (zero-radius OTU).
# ZOTUs are valid OTUs for diversity analysis etc., though the interpretation of
# the results is a bit different from the usual 97% OTUs. For example, it is
# expected that one species may have more than one ZOTU, and with 97% OTUs it is
# expected than an OTU may have more than one species.
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -unoise3 ${Prefix}_derep.fa -zotus ${Prefix}_zOTUs_usearch_unparsed.fa #-relabel OTU #-minampsize 8
# usearch may need "zOUTs be names as OTUs in later steps."
cat ${Prefix}_zOTUs_usearch_unparsed.fa | sed -e 's/Zotu/OTU/' > ${Prefix}_cluster_usearch_zOTUs.fa

#usearch -unoise ${Prefix}.sorted.fa -tabbedout ${Prefix}.txt -faout ${Prefix}.otus.fa -relabel OTU #-minampsize 8

#perl -pi -e 's/uniq.*/OTU . ++$n/ge' ${Prefix}.otus.fa

#rm ${Prefix}.sorted.fa
mv ${Prefix}_derep.fa $CurDir/$OutDir/.
mv ${Prefix}_cluster_vsearch_OTUs.fa $CurDir/$OutDir/.
mv ${Prefix}_cluster_vsearch_zOTUs.fa $CurDir/$OutDir/.
mv ${Prefix}_cluster_usearch_zOTUs.fa $CurDir/$OutDir/.
mv ${Prefix}_cluster_usearch_OTUs.fa $CurDir/$OutDir/.
