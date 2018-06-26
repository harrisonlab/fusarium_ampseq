
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G


QueryReads=$1
RefDb=$2
OtuType=$3
Prefix=$4
OutDir=$5

CurDir=$PWD
WorkDir=${TMPDIR}/quantification
mkdir -p $WorkDir
cd $WorkDir

  # ReadsF=
  # ReadsR=
  # OutF=$(basename ${ReadsF%.fq}_renamed.fa)
  # OutR=$(basename ${ReadsF%.fq}_renamed.fa)
  #
  # #Variable Prefix and SL are given to awk to rename reads and trim them
  # cat $ReadsF | awk -v S="$Prefix" -v SL="$SL" -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){$1=substr($1,(SL+1));print $1}}' > $OutDir/$OutF
  # cat $ReadsR | awk -v S="$Prefix" -v SL="$SL" -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){$1=substr($1,(SL+1));print $1}}' > $OutDir/$OutR
  # # The submit cat file script step of concatenating reads was not performed.
  # # Presumably the step merges filtered and ambiguous reads.

ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -otutab $CurDir/$QueryReads -db $CurDir/$RefDb -strand plus -id 0.97 -biomout ${Prefix}_${OtuType}_table.biom -otutabout ${Prefix}_${OtuType}_table.txt -notmatched ${Prefix}_${OtuType}_nomatch.fa -userout ${Prefix}_${OtuType}_hits.out -userfields query+target
# cat ${Prefix}_${OtuType}_hits.out | cut -f2 | sort | uniq -c | sort -nr > ${Prefix}_${OtuType}_quant.txt
#
# while read Line; do
#   Reads=$(echo $Line | sed -r "s/^\s+//g" | cut -f1 -d ' ')
#   OTU=$(echo $Line | sed -r "s/^\s+//g" | cut -f2 -d ' ')
#   Genus=$(cat $CurDir/$RefTaxa | grep -w "$OTU" | rev | cut -f4 -d ',' | rev )
#   Species=$(cat $CurDir/$RefTaxa | grep -w "$OTU" | rev | cut -f2 -d ',' | rev )
#   printf "$OTU\t$Genus\t$Species\t$Reads\n"
# done < ${Prefix}_${OtuType}_quant.txt > ${Prefix}_${OtuType}_quant.tsv

$ProgDir/usearch -otutab_norm ${Prefix}_${OtuType}_table.txt -sample_size 10000 -output ${Prefix}_${OtuType}_table_norm.txt

mkdir -p $CurDir/$OutDir
cp ${Prefix}_${OtuType}* $CurDir/$OutDir/.
