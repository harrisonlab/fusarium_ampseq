#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

R1=$1
R2=$2
OutDir=$3
# shift
# shift
# PrimerSets=$@

CurDir=$PWD
WorkDir=${TMPDIR}/demulti
mkdir -p $WorkDir
cd $WorkDir

F=$(basename $R1 | cut -f1 -d '.')
R=$(basename $R2 | cut -f1 -d '.')
# F=$(echo $R1|awk -F"/" '{print $NF}')
# R=$(echo $R2|awk -F"/" '{print $NF}')

# OUTDIR=$(echo $R1 |sed "s/$F//" | sed 's/raw_dna/demulti_dna/g')
echo "Output directory:"
echo $CurDir/${OutDir}
mkdir -p $CurDir/"$OutDir"


mkfifo $F.fa
mkfifo $R.fa

zcat -f -- $CurDir/$R1 > $F.fa &
zcat -f -- $CurDir/$R2 > $R.fa &



#ITS
P1=ITS
P1F=GTGAATCATCGAATCTTTGAACGC
P1R=CCGCTTATTGATATGCTTAA[AG]TTCAG
# TEF
P2=TEF
P2F=GGTCACTTGATCTACCAGTGCG
P2R=CCCA[AG]GCGTACTTGAAGRAAC
# SIX13
P3=SIX13
P3F=GCTACTCAAAGTCGTGGACGAG
P3R=GGCAATATATTCCGTCCATTCTTGG
# T2 FOCg17143 OG12981
P4=OG12981
P4F=CACTTCCTCACTTACTTTACCACTCC
P4R=GTCATCGCAATCGCCKTCCG
# T4 orthogroup 13890
P5=OG13890
P5F=GCTGTCTTATCACTTATCAGCCTTG
P5R=CGGTCTGATTTGGTGTCCAGTCG

# ProgDir=/home/armita/git_repos/emr_repos/scripts/Metabarcoding_pipeline/scripts
# # $ProgDir/demulti_v2.pl tmpF.fq tmpR.fq $P1F $P1R $P2F $P2R $P3F $P3R $P4F $P4R $P5F $P5R
# $ProgDir/demulti_v2.pl $F.fa $R.fa $P1F $P1R $P2F $P2R $P3F $P3R $P4F $P4R $P5F $P5R

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
# $ProgDir/demulti.py --FastqF tmpF.fq --FastqR tmpR.fq --primer_loci $P1 $P2 $P3 $P4 $P5 --primersF $P1F $P2F $P3F $P4F $P5F --primersR $P1R $P2R $P3R $P4R $P5R
$ProgDir/demulti.py --FastqF $F.fa --FastqR $R.fa --primer_loci $P1 $P2 $P3 $P4 $P5 --primersF $P1F $P2F $P3F $P4F $P5F --primersR $P1R $P2R $P3R $P4R $P5R  2>&1 | tee ${Prefix}_demulti.log

rm $F.fa $R.fa

cp *.fq $CurDir/$OutDir/.
cp *.fastq $CurDir/$OutDir/.
