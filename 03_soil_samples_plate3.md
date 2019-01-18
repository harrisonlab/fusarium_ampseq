# fusarium_ampseq
Commands used in analysis of amplicon sequence data of Fusarium infested soil samples and artificial mixes.

All of this work was done in the directory:
```bash
  cd /home/groups/harrisonlab/project_files/fusarium_ampseq
```

# Data transfer

Samples sequenced were:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/181207_M04465_0090_000000000-BV5PW/Data/Intensities/BaseCalls
  ls $RawDatDir
```

Sample names referred to the following:

Loci:
16S, ITS, TEF, T1 (SIX13), T4 (OG13890), T6 (), T8, with combined 16S and ITS samples run using a single barcode and TEF using a single barcode with multiplexed T1, T4, T6  or with mtuliplexed T1, T4, T6 and T8.


## Sample assignment

Details of sample ID were listed and copied into an excel file. Run information was added to the file and it was saved as a TSV file. The contents of this file were opened using a TextEdit and then copied and pasted into a nano-generated text file on the cluster.


```bash
RawDatDir=/data/seq_data/miseq/2018/RAW/181207_M04465_0090_000000000-BV5PW/Data/Intensities/BaseCalls
for File in $(ls $RawDatDir/*.fastq.gz); do
  basename $File
done

cd /home/groups/harrisonlab/project_files/fusarium_ampseq
mkdir data/plate3

nano data/plate3/plate3_sample_assignment.txt
```

Raw sequencing data was symbolicly linked to the working directory:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/181207_M04465_0090_000000000-BV5PW/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
  AssignmentFile=$(ls data/plate3/plate3_sample_assignment.txt)
  for ReadF in $(ls $RawDatDir/*_L001_R1_001.fastq.gz | grep -v 'Undetermined'); do
    ReadR=$(echo $ReadF | sed 's/_L001_R1_001.fastq.gz/_L001_R2_001.fastq.gz/g')

    Name=$(basename $ReadF)
    Locus=$(cat $AssignmentFile | grep "$Name" | cut -f2)
    Field=$(cat $AssignmentFile | grep "$Name" | cut -f3)
    Block=$(cat $AssignmentFile | grep "$Name" | cut -f4)
    TechRep=$(cat $AssignmentFile | grep "$Name" | cut -f5)
    echo "$Locus - $Field - $Block - $TechRep"
    # OutDir=/data2/scratch2/armita/fusarium_ampseq/raw_dna/plate3/paired/$Locus/$Pool/$Dilution/$TechRep
    OutDir=$PWD/raw_dna/plate3/paired/$Locus/$Field/$Block/$TechRep
    mkdir -p $OutDir/F
    # cp $ReadF $OutDir/F/.
    # mkdir -p $OutDir/R
    # cp $ReadR $OutDir/R/.
    cp -s $ReadF $OutDir/F/.
    mkdir -p $OutDir/R
    cp -s $ReadR $OutDir/R/.
  done
```

# Metabarcoding analysis

These commands are based upon the workflow described at:
https://github.com/eastmallingresearch/Metabarcoding_pipeline

### Set pipeline variables

```bash
# set MBPL variable to pipeline folder
MBPL=/home/deakig/metabarcoding_pipeline
# to set permanetly for future (bash) shell sessions (be careful with this, if you have settings in ~/.profile they will no longer load)
# echo export MBPL=~/metabarcoding_pipeline >>~/.bash_profile
```

Prerequesits:
*HMM Preperation for ITS analysis*
Greg has prepared hmm models to identify the 5.8s, ssu and lsu regions of
the ITS region. These can be removed from the sequences before blast Searching
for isolate taxonomy.
*Taxonomy reference databases*
Greg has downloaded the Unite V7 (fungi) and RDP trainset 15 (bacteria) reference databases from
http://drive5.com/usearch/manual/utax_downloads.html
Configuration has been tested with usearch8.1 (probably works the same for 10)


Greg has also built an oomycota database combining three sets of data; 1) a subset (stamenopiles) of the silva_ssu database https://www.arb-silva.de/browser/, 2) a supplied 3rd party database 3) and a usearch specific Unite database (https://unite.ut.ee/sh_files/utax_reference_dataset_28.06.2017.zip)

NOTE: It is now possible to download just the stramenopiles subset from silva

# Analysis of individual loci

## Demultiplexing

This script demultiplexs mixed (e.g. ITS and 16S) libraries based on the primer sequence. Any sequence which has mismatches is written to ambiguous.fq (f & r seperately). Primer sequences
are detailed in the submission wrapper.
*Note* Regex are used to describe degenerate bases in the primer.

Run below to demultiplex:

```bash
for DataDir in $(ls -d raw_dna/plate3/paired/*/*/*/*); do
Jobs=$(qstat | grep 'submit_dem' | grep 'r' | grep 'qw' | wc -l)
while [ $Jobs -gt 5 ]; do
sleep 10s
printf "."
Jobs=$(qstat | grep 'submit_dem' | grep 'r' | grep 'qw' | wc -l)
done
printf "\n"
R1=$(ls $DataDir/F/*.fastq.gz)
R2=$(ls $DataDir/R/*.fastq.gz)
echo $DataDir
echo $(basename $R1)
echo $(basename $R2)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
OutDir=demulti_dna/plate3/$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
qsub $ProgDir/submit_demulti.sh $R1 $R2 ${OutDir}
done
```


 Summise reads demultiplexed:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\t16S\tITS\tTef\tSIX13\tT4\tT6\tT8\tAmbiguous\n" > demulti_dna/plate3/demultiplex_summary.tsv
for RunDir in $(ls -d demulti_dna/plate3/*/*/*/*); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  T0=$(ls $RunDir/*_R1_*16S.fq)
  T0Lines=$(cat $T0 | awk '{s++}END{print s/4}')
  ITS=$(ls $RunDir/*_R1_*ITS.fq)
  ItsLines=$(cat $ITS | awk '{s++}END{print s/4}')
  TEF=$(ls $RunDir/*_R1_*TEF.fq)
  TefLines=$(cat $TEF | awk '{s++}END{print s/4}')
  Six13=$(ls $RunDir/*_R1_*SIX13.fq)
  Six13Lines=$(cat $Six13 | wc -l)
  # T2=$(ls $RunDir/*OG12981.fq)
  # T2Lines=$(cat $T2 | wc -l)
  T4=$(ls $RunDir/*_R1_*OG13890.fq)
  T4Lines=$(cat $T4 | awk '{s++}END{print s/4}')
  T6=$(ls $RunDir/*_R1_*OG4952.fq)
  T6Lines=$(cat $T6 | awk '{s++}END{print s/4}')
  T8=$(ls $RunDir/*_R1_*OG13397.fq)
  T8Lines=$(cat $T8 | awk '{s++}END{print s/4}')
  Ambiguous=$(ls $RunDir/*_R1_*ambiguous.fq)
  AmbLines=$(cat $Ambiguous | awk '{s++}END{print s/4}')
  printf "$RunName\t$Locus\t$Pool\t$Dilution\t$Rep\t$T0Lines\t$ItsLines\t$TefLines\t$Six13Lines\t$T4Lines\t$T6Lines\t$T8Lines\t$AmbLines\n"
done >> demulti_dna/plate3/demultiplex_summary.tsv
ls $PWD/demulti_dna/plate3/demultiplex_summary.tsv
```

From this data thresholding values of cross-contamination as a result of illumina
adapter read hopping were determined. For each row in the dataset, the most
abundant locus was assumed to be the target locus and reads attributed to other
loci assumed to be contaminant reads. Reads from another experiment with the same
locus could be contaminating the sample. The 2nd most abundant locus was
identified for each run (the highest contaminant locus) and the maximum value
identified across the entire plate. A threshold for a minimum abundance of reads
attributed to an OTU was set at this value.

In the case of our plate, this was:
```bash
Threshold=
```

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.


```bash
for DataDir in $(ls -d demulti_dna/plate3/mix-A/*/*/*); do
  for R1 in $(ls $DataDir/*R1*.fq | grep -v 'ambiguous' | grep -e '16S.fq' -e 'ITS.fq'); do
    R2=$(echo $R1 | sed 's/_R1/_R2/g')
    echo $(basename $R1)
    echo $(basename $R2)
    Locus=$(echo $R1 | rev | cut -f1 -d '_' | rev | tr -d '.fq')
    Prefix=$(echo $DataDir | cut -f3,4,5,6 -d '/' | sed 's&/&_&g')"_${Locus}"
    OutDir="processed_dna/"$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
    Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    while [ $Jobs -gt 5 ]; do
    sleep 10s
    printf "."
    Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    done
    printf "\n"
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
  done
done
for DataDir in $(ls -d demulti_dna/plate3/mix-B/*/*/*); do
  for R1 in $(ls $DataDir/*R1*.fq | grep -v 'ambiguous' | grep -e 'TEF.fq' -e 'SIX13.fq' -e 'OG13890.fq' -e 'OG4952.fq'); do
    R2=$(echo $R1 | sed 's/_R1/_R2/g')
    echo $(basename $R1)
    echo $(basename $R2)
    Locus=$(echo $R1 | rev | cut -f1 -d '_' | rev | tr -d '.fq')
    Prefix=$(echo $DataDir | cut -f3,4,5,6 -d '/' | sed 's&/&_&g')"_${Locus}"
    OutDir="processed_dna/"$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
    Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    while [ $Jobs -gt 5 ]; do
    sleep 10s
    printf "."
    Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    done
    printf "\n"
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
  done
done
for DataDir in $(ls -d demulti_dna/plate3/mix-C/*/*/*); do
  for R1 in $(ls $DataDir/*R1*.fq | grep -v 'ambiguous' | grep -e 'TEF.fq' -e 'SIX13.fq' -e 'OG13890.fq' -e 'OG4952.fq' -e 'OG13397.fq'); do
    R2=$(echo $R1 | sed 's/_R1/_R2/g')
    echo $(basename $R1)
    echo $(basename $R2)
    Locus=$(echo $R1 | rev | cut -f1 -d '_' | rev | tr -d '.fq')
    Prefix=$(echo $DataDir | cut -f3,4,5,6 -d '/' | sed 's&/&_&g')"_${Locus}"
    OutDir="processed_dna/"$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
    Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    while [ $Jobs -gt 5 ]; do
    sleep 10s
    printf "."
    Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    done
    printf "\n"
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
  done
done
```


 Summise reads merged:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\t16S\tITS\tTef\tSIX13\tT4\tT6\tT8\tAmbiguous\n" > processed_dna/plate3/merged_summary2.tsv
for RunDir in $(ls -d demulti_dna/plate3/*/*/*/*); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName="processed_dna/plate3/$Locus/$Pool/$Dilution/merged/${Locus}_${Pool}_${Dilution}_${Rep}*"
  # RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  T0Lines=$(cat ${RunName}*_*16S_prefilter.fa | grep '>'| wc -l)
  ItsLines=$(cat ${RunName}*_*ITS_prefilter.fa | grep '>'| wc -l)
  TefLines=$(cat ${RunName}*_*TEF_prefilter.fa | grep '>'| wc -l)
  Six13=$(ls $RunDir/*_R1_*SIX13.fq)
  Six13Lines=$(cat ${RunName}*_*SIX13_prefilter.fa | grep '>'| wc -l)
  T4Lines=$(cat ${RunName}*_*OG13890_prefilter.fa | grep '>'| wc -l)
  T6Lines=$(cat ${RunName}*_*OG4952_prefilter.fa | grep '>'| wc -l)
  T8Lines=$(cat ${RunName}*_*OG13397_prefilter.fa | grep '>'| wc -l)
  printf "${Locus}_${Pool}_${Dilution}_${Rep}\t$Locus\t$Pool\t$Dilution\t$Rep\t$T0Lines\t$ItsLines\t$TefLines\t$Six13Lines\t$T4Lines\t$T6Lines\t$T8Lines\n"
done >> processed_dna/plate3/merged_summary2.tsv
ls $PWD/processed_dna/plate3/merged_summary2.tsv
```


Summise reads filtered:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\t16S\tITS\tTef\tSIX13\tT4\tT6\tT8\tAmbiguous\n" > processed_dna/plate3/filtered_summary2.tsv
for RunDir in $(ls -d demulti_dna/plate3/*/*/*/*); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName="processed_dna/plate3/$Locus/$Pool/$Dilution/filtered/${Locus}_${Pool}_${Dilution}_${Rep}*"
  # RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  T0Lines=$(cat ${RunName}*_*16S.filtered.fa | grep '>'| wc -l)
  ItsLines=$(cat ${RunName}*_*ITS.filtered.fa | grep '>'| wc -l)
  TefLines=$(cat ${RunName}*_*TEF.filtered.fa | grep '>'| wc -l)
  Six13Lines=$(cat ${RunName}*_*SIX13.filtered.fa | grep '>'| wc -l)
  T4Lines=$(cat ${RunName}*_*OG13890.filtered.fa | grep '>'| wc -l)
  T6Lines=$(cat ${RunName}*_*OG4952.filtered.fa | grep '>'| wc -l)
  T8Lines=$(cat ${RunName}*_*OG13397.filtered.fa | grep '>'| wc -l)
  printf "${Locus}_${Pool}_${Dilution}_${Rep}\t$Locus\t$Pool\t$Dilution\t$Rep\t$T0Lines\t$ItsLines\t$TefLines\t$Six13Lines\t$T4Lines\t$T6Lines\t$T8Lines\n"
done >> processed_dna/plate3/filtered_summary2.tsv
ls $PWD/processed_dna/plate3/filtered_summary2.tsv
```


## OTU assignment
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting

 * All read files associated with the each locus in the project are concatenated
 * Clustering run on all the data for this locus within the project
 * Quantification can then be performed for each treatment against the total set

```bash
  # Concatenate files
  for Locus in 16S ITS TEF SIX13 OG13397 OG13890 OG4952; do
    OutDir=clustering/plate3/$Locus
    mkdir -p $OutDir
    cat processed_dna/plate3/*/*/*/filtered/*_${Locus}.filtered.fa >> $OutDir/${Locus}_concatenated.temp.fa
    ls -lh $OutDir/${Locus}_concatenated.temp.fa
    # Submit clustering
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/submit_clustering.sh $OutDir/${Locus}_concatenated.temp.fa $OutDir $Locus
  done
```


### Assign taxonomy

https://www.drive5.com/usearch/manual/utax_or_sintax.html
Taxonomy can be assigned using SINTAX or UTAX methods:

```
UTAX and SINTAX have different strengths and weaknesses
SINTAX is brand new so I don't have much experience with it yet (this was written just after version 9 was released). On short 16S tags like V4, SINTAX and RDP have very similar performance. On longer 16S sequences and on ITS sequences, SINTAX is better than RDP. SINTAX is simpler because it doesn't need training, while training UTAX or RDP is quite challenging if you want to use your own database. UTAX is the only algorithm which tries to account for sparse reference data and has the lowest over-classification rate of any algorithm (except possibly the k-nearest-neighbor method in mothur, but knn has low sensitivity in general). However, UTAX sometimes has lower sensitivity than SINTAX to known taxa. Neither algorithm is a clear winner over the other.
```

* Greg uses Sintax as it is A) recomended by the author B) easier to make custom databases with
* An alternative to database-based taxonomy assignment is to used BLAST. It may be worth investigating automated BLAST lookups vs NCBI.

usearch hosts databases for some loci. Silva also has a database for the stremenophiles (oomycetes) that can be downloaded. Otherwise, personal databases must be made.
https://unite.ut.ee/repository.php
https://www.arb-silva.de/browser/

The rdp bacterial 16S database was downloaded:

```bash
qlogin

WorkDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd $WorkDir
OutDir=databases/16S
mkdir -p $OutDir
cd $OutDir

# wget -N https://drive5.com/sintax/rdp_16s_v16.fa.gz
# wget didnt sucessfully download the file so i downloaded the file to my local
# computer from: https://drive5.com/sintax/
# Then copied the file to the cluster:
# scp /Users/armita/Downloads/rdp_16s_v16.fa.gz cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/databases/16S/.

gunzip *.gz
mv rdp_16s_v16.fa rdp_bacterial_16S.fasta
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax rdp_bacterial_16S.fasta -output rdp_bacterial_16S_sintax.udp
cd $WorkDir

logout
```



```bash
for Locus in 16S ITS TEF SIX13 OG13397 OG13890 OG4952; do
  for Type in zOTUs; do
  # for Type in OTUs zOTUs; do
  OtuFa=$(ls clustering/plate3/$Locus/${Locus}_${Type}.fa)
  RefDb=$(ls databases/$Locus/*.udp)
  Prefix=$(basename ${OtuFa%.fa})
  OutDir=$(dirname $OtuFa)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
  # echo "$OtuFa $RefDb $Prefix $OutDir"
  qsub $ProgDir/submit_taxonomy.sh $OtuFa $RefDb $Prefix $OutDir
  done
done
```

## Quantification

```bash
  for Locus in mix-A; do
    for Field in onion daffodil stocks; do
      OutDir=quantified/plate3/$Locus/$Field
      mkdir -p $OutDir
      Threshold=114
      Identity=0.97
      for OtuType in zOTUs; do
        for RefDb in $(ls clustering/plate3/*/*_${OtuType}_taxa.fa | grep -e '16S' -e 'ITS' | grep '16S'); do
          QueryReads=$OutDir/${Locus}_$(basename $RefDb | cut -f1 -d '_')_reads_appended.fa
          cat processed_dna/plate3/$Locus/$Field/*/merged/*$(basename $RefDb | cut -f1 -d '_')_prefilter.fa | cut -f1 -d '.' | sed 's/mix-/mix_/g' > $QueryReads
          Prefix="${Locus}_$(basename $RefDb | cut -f1 -d '_')_$Field_${OtuType}"
          ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
          qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
        done
      done
    done
  done
```


```bash
  for Locus in mix-B; do
    for Field in onion daffodil stocks; do
      OutDir=quantified/plate3/$Locus/$Field
      mkdir -p $OutDir
      Threshold=114
      Identity=1
      for OtuType in zOTUs; do
        # for RefDb in $(ls clustering/plate3/*/*_${OtuType}_taxa.fa | grep -e 'TEF' -e 'SIX13' -e 'OG13890' -e 'OG4952'); do
        for RefDb in $(ls databases/*/*_amplicon_nr_db.fasta | grep -e 'TEF' -e 'SIX13' -e 'OG13890' -e 'OG4952'); do
          Locus2=$(basename $RefDb | cut -f1 -d '_')
          QueryReads=$OutDir/${Locus}_$(basename $RefDb | cut -f1 -d '_')_reads_appended.fa
          # cat processed_dna/plate3/$Locus/$Field/*/merged/*$(basename $RefDb | cut -f1 -d '_')_prefilter.fa | cut -f1 -d '.' | sed 's/mix-/mix_/g' > $QueryReads
          cat processed_dna/plate3/$Locus/$Field/*/filtered/*$(basename $RefDb | cut -f1 -d '_').filtered.fa | cut -f1 -d '.' | sed "s/${Locus2}.*/${Locus2}/g" | sed 's/mix-/mix_/g' > $QueryReads
          Prefix="${Locus}_$(basename $RefDb | cut -f1 -d '_')_$Field_${OtuType}"
          ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
          qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
        done
      done
    done
  done
```

```bash
  for Locus in mix-C; do
    for Field in onion daffodil stocks; do
      OutDir=quantified/plate3/$Locus/$Field
      mkdir -p $OutDir
      Threshold=114
      Identity=1
      for OtuType in zOTUs; do
        for RefDb in $(ls databases/*/*_amplicon_nr_db.fasta | grep -e 'TEF' -e 'SIX13' -e 'OG13890' -e 'OG4952' -e 'OG13397'); do
          Locus2=$(basename $RefDb | cut -f1 -d '_')
          QueryReads=$OutDir/${Locus}_$(basename $RefDb | cut -f1 -d '_')_reads_appended.fa
          cat processed_dna/plate3/$Locus/$Field/*/filtered/*$(basename $RefDb | cut -f1 -d '_').filtered.fa | cut -f1 -d '.' | sed "s/${Locus2}.*/${Locus2}/g" | sed 's/mix-/mix_/g' > $QueryReads
          Prefix="${Locus}_$(basename $RefDb | cut -f1 -d '_')_$Field_${OtuType}"
          ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
          echo "qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity"
        done
      done
    done
  done
```



```bash
  rm quantified/plate3/*/*/*.fa
  rm quantified/plate3/*/*/*_hits.out
  ls -d $PWD/quantified/plate3
```



```bash
cd /Users/armita/Downloads/AHDB_new/plate3

for File in $(ls quantified/mix-*/*/mix-*_TEF_zOTUs_zOTUs_table.txt | grep -v 'daffodil'); do
Prefix=$(basename $File | sed 's/_zOTUs_zOTUs_table.txt//g')
OutDir=$(dirname $File | sed 's/quantified/plots/g')
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_TEF.r --OTU_table $File --prefix $OutDir/$Prefix
done

for File in $(ls quantified/mix-*/*/mix-*_*_zOTUs_zOTUs_table.txt | grep -v 'TEF' | grep -v 'daffodil' | grep -e 'OG4952' -e 'OG13397'); do
Prefix=$(basename $File | sed 's/_zOTUs_zOTUs_table.txt//g')
OutDir=$(dirname $File | sed 's/quantified/plots/g')
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_SIX13.r --OTU_table $File --prefix $OutDir/$Prefix
done
```

<!-- ```bash
for File in $(ls quantified/mix-C/*/mix-C_TEF_zOTUs_zOTUs_table.txt); do
Prefix=$(basename $File | sed 's/_zOTUs_zOTUs_table.txt//g')
OutDir=$(dirname $File | sed 's/quantified/plots/g')
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_TEF.r --OTU_table $File --prefix $OutDir/$Prefix
done

for File in $(ls quantified/mix-C/*/mix-C_*_zOTUs_zOTUs_table.txt | grep -v 'TEF'); do
Prefix=$(basename $File | sed 's/_zOTUs_zOTUs_table.txt//g')
OutDir=$(dirname $File | sed 's/quantified/plots/g')
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_SIX13.r --OTU_table $File --prefix $OutDir/$Prefix
done
```
 -->

```bash
for File in $(ls quantified/mix-A/*/mix-A_ITS_zOTUs_zOTUs_table.txt | grep -v 'daffodil'); do
Prefix=$(basename $File | sed 's/_zOTUs_zOTUs_table.txt//g')
OutDir=$(dirname $File | sed 's/quantified/plots/g')
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_ITS.r --OTU_table $File --prefix $OutDir/$Prefix
done

for File in $(ls quantified/mix-A/*/mix-A_16S_zOTUs_zOTUs_table.txt | grep -v 'daffodil'); do
Prefix=$(basename $File | sed 's/_zOTUs_zOTUs_table.txt//g')
OutDir=$(dirname $File | sed 's/quantified/plots/g')
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_16S.r --OTU_table $File --prefix $OutDir/$Prefix
done
```

```bash
cd /Users/armita/Downloads/AHDB_new/plate3

OnionTab=$(ls quantified/mix-A/onion/mix-A_16S_zOTUs_zOTUs_table.txt)
DaffodilTab=$(ls quantified/mix-A/daffodil/mix-A_16S_zOTUs_zOTUs_table.txt)
StocksTab=$(ls quantified/mix-A/stocks/mix-A_16S_zOTUs_zOTUs_table.txt)
Prefix="mix-A_16S_by_field"
OutDir="plots/mix-A"
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_16s_all_fields.r --OTU_table_onion $OnionTab --OTU_table_daffodil $DaffodilTab --OTU_table_stocks $StocksTab --prefix $OutDir/$Prefix

OnionTab=$(ls quantified/mix-A/onion/mix-A_ITS_zOTUs_zOTUs_table.txt)
DaffodilTab=$(ls quantified/mix-A/daffodil/mix-A_ITS_zOTUs_zOTUs_table.txt)
StocksTab=$(ls quantified/mix-A/stocks/mix-A_ITS_zOTUs_zOTUs_table.txt)
Prefix="mix-A_ITS_by_field"
OutDir="plots/mix-A"
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_16s_all_fields.r --OTU_table_onion $OnionTab --OTU_table_daffodil $DaffodilTab --OTU_table_stocks $StocksTab --prefix $OutDir/$Prefix

OnionTab=$(ls quantified/mix-B/onion/mix-B_TEF_zOTUs_zOTUs_table.txt)
DaffodilTab=$(ls quantified/mix-B/daffodil/mix-B_TEF_zOTUs_zOTUs_table.txt)
StocksTab=$(ls quantified/mix-B/stocks/mix-B_TEF_zOTUs_zOTUs_table.txt)
Prefix="mix-B_TEF_by_field"
OutDir="plots/mix-B"
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_TEF_all_fields.r --OTU_table_onion $OnionTab --OTU_table_daffodil $DaffodilTab --OTU_table_stocks $StocksTab --prefix $OutDir/$Prefix

for Locus in OG4952; do
OnionTab=$(ls quantified/mix-B/onion/mix-B_${Locus}_zOTUs_zOTUs_table.txt)
DaffodilTab=$(ls quantified/mix-B/daffodil/mix-B_${Locus}_zOTUs_zOTUs_table.txt)
StocksTab=$(ls quantified/mix-B/stocks/mix-B_${Locus}_zOTUs_zOTUs_table.txt)
Prefix="mix-B_${Locus}_by_field"
OutDir="plots/mix-B"
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_fsp_loci_all_fields.r --OTU_table_onion $OnionTab --OTU_table_daffodil $DaffodilTab --OTU_table_stocks $StocksTab --prefix $OutDir/$Prefix
done

for Locus in SIX13 OG13890; do
OnionTab=$(ls quantified/mix-B/onion/mix-B_${Locus}_zOTUs_zOTUs_table.txt)
DaffodilTab=$(ls quantified/mix-B/daffodil/mix-B_${Locus}_zOTUs_zOTUs_table.txt)
StocksTab=$(ls quantified/mix-B/stocks/mix-B_${Locus}_zOTUs_zOTUs_table.txt)
Prefix="mix-B_${Locus}_by_field_total"
OutDir="plots/mix-B"
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_fsp_loci_all_fields_total_counts.r --OTU_table_onion $OnionTab --OTU_table_daffodil $DaffodilTab --OTU_table_stocks $StocksTab --prefix $OutDir/$Prefix
done

for Locus in OG13397; do
OnionTab=$(ls quantified/mix-C/onion/mix-C_${Locus}_zOTUs_zOTUs_table.txt)
DaffodilTab=$(ls quantified/mix-C/daffodil/mix-C_${Locus}_zOTUs_zOTUs_table.txt)
StocksTab=$(ls quantified/mix-C/stocks/mix-C_${Locus}_zOTUs_zOTUs_table.txt)
Prefix="mix-C_${Locus}_by_field"
OutDir="plots/mix-C"
mkdir -p $OutDir
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/plate3
$ProgDir/plot_plate3_fsp_loci_all_fields.r --OTU_table_onion $OnionTab --OTU_table_daffodil $DaffodilTab --OTU_table_stocks $StocksTab --prefix $OutDir/$Prefix
done
```
