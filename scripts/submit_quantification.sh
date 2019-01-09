
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=2G
#$ -l h=blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace


QueryReads=$1
RefDb=$2
OtuType=$3
Prefix=$4
OutDir=$5
Threshold=$6
Identity=$7

CurDir=$PWD
WorkDir=${TMPDIR}/quantification
mkdir -p $WorkDir
cd $WorkDir

# quantify taxa
ProgDir=/home/deakig/usr/local/bin
# $ProgDir/usearch -otutab $CurDir/$QueryReads -db $CurDir/$RefDb -strand plus -id 0.97 -biomout ${Prefix}_${OtuType}_table.biom -otutabout ${Prefix}_${OtuType}_table.txt -notmatched ${Prefix}_${OtuType}_nomatch.fa -userout ${Prefix}_${OtuType}_hits.out -userfields query+target
$ProgDir/usearch -otutab $CurDir/$QueryReads -db $CurDir/$RefDb -strand plus -id $Identity -biomout ${Prefix}_${OtuType}_table.biom -otutabout ${Prefix}_${OtuType}_table.txt -notmatched ${Prefix}_${OtuType}_nomatch.fa -userout ${Prefix}_${OtuType}_hits.out -userfields query+target
# Normalise results to 10000 reads
$ProgDir/usearch -otutab_norm ${Prefix}_${OtuType}_table.txt -sample_size 10000 -output ${Prefix}_${OtuType}_table_norm.txt

# Combine reads by species and filter taxa with reads fewer than a given threshold
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
$ProgDir/filter_OTU_tables.py --table ${Prefix}_${OtuType}_table.txt --threshold 0 --prefix ${Prefix}_${OtuType}
$ProgDir/filter_OTU_tables.py --table ${Prefix}_${OtuType}_table.txt --threshold $Threshold --prefix ${Prefix}_${OtuType}_thresholded
# Normalise these filtered samples
cat ${Prefix}_${OtuType}_thresholded_table_by_spp.txt | sed "s/^Species/#OTU ID/g" > tmp.txt
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -otutab_norm tmp.txt -sample_size 10000 -output tmp2.txt
cat tmp2.txt | sed "s/^#OTU ID/Species/g"  > ${Prefix}_${OtuType}_thresholded_norm_table_by_spp.txt

cat ${Prefix}_${OtuType}_thresholded_table_by_genus.txt | sed "s/^Genus/#OTU ID/g" > tmp.txt
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -otutab_norm tmp.txt -sample_size 10000 -output tmp2.txt
cat tmp2.txt | sed "s/^#OTU ID/Genus/g"  > ${Prefix}_${OtuType}_thresholded_norm_table_by_genus.txt

mkdir -p $CurDir/$OutDir
cp ${Prefix}_${OtuType}* $CurDir/$OutDir/.
