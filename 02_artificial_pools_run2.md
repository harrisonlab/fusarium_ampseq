# fusarium_ampseq
Commands used in analysis of amplicon sequence data of Fusarium infested soil samples and artificial mixes.

All of this work was done in the directory:
```bash
  cd /home/groups/harrisonlab/project_files/fusarium_ampseq
```

# Data transfer

Samples sequenced were:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/181126_M04465_0089_000000000-BV5TT/Data/Intensities/BaseCalls
  ls $RawDatDir
```

Sample names referred to the following:

Loci:
ITS, TEF, T1 (SIX13), T2 (FoC_g17143 / OG12981), T4 (OG13890), PCR1 (pooling of isolate pool 1 and 2 with all loci at equimolar concentrations, pooled for multiplex PCR of all loci), PCR2 (pooling of isolate pool 1 with all loci at equimolar concentrations, pooled following individual PCRs of all loci before adapter ligation PCR)

Pool:
Samples 1-24 were performed on isolate pool 2, samples 25-28 were perfomed on all isolates from both pools, samples 29-80 were performed on isolate pool 1.

Dilution:
a = equimolar of all isolates
b = 1 in 10 dilution of 1 or 2 isolates
c = 1 in 100 dilution of 1 or 2 isolates

Sample assignment:
```bash
cd /home/groups/harrisonlab/project_files/fusarium_ampseq
mkdir data/plate2c
# run from local machine
# cat /Users/armita/Downloads/Plate2c_sample_assignment.txt | sed 's/rep1/rep1\n/g'| sed 's/rep2/rep2\n/g' | sed 's/rep3/rep3\n/g'| sed -e "s/\tpool/\tpool\n/g" | sed 's/\/pool)/pool)\n/g' | tr -d '\r' > /Users/armita/Downloads/Plate2c_sample_assignment_ed.txt

# scp /Users/armita/Downloads/Plate2c_sample_assignment.txt cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/data/plate2c/Plate2c_sample_assignment2.txt
# cat data/plate2c/Plate2c_sample_assignment2.txt | sed 's/rep1/rep1\n/g'| sed 's/rep2/rep2\n/g' | sed 's/rep3/rep3\n/g'| sed -e "s/\tpool/\tpool\n/g" | sed "s/\/pool)/pool)\n/g" | tr -d '\r' > data/plate2c/Plate2c_sample_assignment_ed.txt

# The tsv file on my local computer was opened in notepad and copied and pasted
# into a document on the cluster. This removed window newline characters.
nano data/plate2c/Plate2c_sample_assignment.txt

```

Raw sequencing data was symbolicly linked to the working directory:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/181126_M04465_0089_000000000-BV5TT/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
  AssignmentFile=$(ls data/plate2c/Plate2c_sample_assignment.txt)
  for ReadF in $(ls $RawDatDir/*_L001_R1_001.fastq.gz | grep -v 'Undetermined'); do
    ReadR=$(echo $ReadF | sed 's/_L001_R1_001.fastq.gz/_L001_R2_001.fastq.gz/g')

    Name=$(basename $ReadF)
    Locus=$(cat $AssignmentFile | grep "$Name" | cut -f2)
    Pool=$(cat $AssignmentFile | grep "$Name" | cut -f3)
    Dilution=$(cat $AssignmentFile | grep "$Name" | cut -f4)
    TechRep=$(cat $AssignmentFile | grep "$Name" | cut -f5)
    echo "$Locus - $Pool - $Dilution - $TechRep"
    # OutDir=/data2/scratch2/armita/fusarium_ampseq/raw_dna/plate2c/paired/$Locus/$Pool/$Dilution/$TechRep
    OutDir=$PWD/raw_dna/plate2c/paired/$Locus/$Pool/$Dilution/$TechRep
    mkdir -p $OutDir/F
    cp $ReadF $OutDir/F/.
    mkdir -p $OutDir/R
    cp $ReadR $OutDir/R/.
    # cp -s $ReadF $OutDir/F/.
    # mkdir -p $OutDir/R
    # cp -s $ReadR $OutDir/R/.
  done
```


<!--
#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
  # cd /data2/scratch2/armita/fusarium_ampseq
  for RawData in $(ls raw_dna/plate2c/paired/*/*/*/*/*/*.fastq.gz); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    Jobs=$(qstat | grep 'run_fastq' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 1m
    printf "."
    Jobs=$(qstat | grep 'run_fastq' | grep 'qw' | wc -l)
    done		
    printf "\n"
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```

```bash
cd /data2/scratch2/armita/fusarium_ampseq
for StrainPath in $(ls -d raw_dna/plate2c/paired/*/*/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
done		
printf "\n"
ReadsF=$(ls $StrainPath/F/*.fastq*)
ReadsR=$(ls $StrainPath/R/*.fastq*)
echo $ReadsF
echo $ReadsR
OutDir=$(echo $StrainPath | sed 's/raw_dna/qc_dna/g')
qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA $OutDir
done
```

The number of reads for each treatment was determined:

```bash
cd /data2/scratch2/armita/fusarium_ampseq
for F_read in $(ls qc_dna/plate2c/paired/*/*/*/*/F/*fq.gz); do
  Locus=$(echo $F_read | cut -f3 -d '/')
  Pool=$(echo $F_read | cut -f4 -d '/')
  Dilution=$(echo $F_read | cut -f5 -d '/')
  TechRep=$(echo $F_read | cut -f6 -d '/')
  ReadCount=$(cat $F_read | gunzip -cf | awk '{s++}END{print s/4}')
  printf "$Locus\t$Pool\t$Dilution\t$TechRep\t$ReadCount\n"
done > qc_dna/reads_per_sample.tsv
``` -->


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
# cd /data2/scratch2/armita/fusarium_ampseq
for DataDir in $(ls -d raw_dna/plate2c/paired/*/*/*/*); do
# for DataDir in $(ls -d raw_dna/plate2c/paired/*/*/*/* | grep -e 'ITS/soil_pathogens/dilution_100x/rep3' -e 'ITS/soil_pathogens/dilution_10x/pool' -e 'og13397/Fusarium_spp/dilution_100x/rep1'); do
Jobs=$(qstat | grep 'submit_dem' | grep 'r' | grep 'qw' | wc -l)
while [ $Jobs -gt 5 ]; do
sleep 10s
printf "."
Jobs=$(qstat | grep 'submit_dem' | grep 'r' | grep 'qw' | wc -l)
done
printf "\n"
# WorkDir=/data2/scratch2/armita/fusarium_ampseq
R1=$(ls $DataDir/F/*.fastq.gz)
R2=$(ls $DataDir/R/*.fastq.gz)
echo $DataDir
echo $(basename $R1)
echo $(basename $R2)
# ls -lh $R1
# ls -lh $R2
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
OutDir=demulti_dna/plate2c/$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
qsub $ProgDir/submit_demulti.sh $R1 $R2 ${OutDir}
done
```

Summise reads demultiplexed:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\tITS\tTef\tSIX13\tT4\tT6\tT8\tAmbiguous\n" > demulti_dna/plate2c/demultiplex_summary.tsv
for RunDir in $(ls -d demulti_dna/plate2c/*/*/*/*); do
 Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
 Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
 Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
 Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
 RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
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
 printf "$RunName\t$Locus\t$Pool\t$Dilution\t$Rep\t$ItsLines\t$TefLines\t$Six13Lines\t$T4Lines\t$T6Lines\t$T8Lines\t$AmbLines\n"
done >> demulti_dna/plate2c/demultiplex_summary.tsv
ls $PWD/demulti_dna/plate2c/demultiplex_summary.tsv
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
Threshold=210
```


<!-- ### Ambiguous data
Ambiguous data should not be used for OTU clustering/denoising, but it can be counted in the OTU tables.
Would require converting to FASTA with approprite labels  - the below should do this will join PE and remove adapters/primers
```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c AMBIGpre \
 "$PROJECT_FOLDER/data/$RUN/ambiguous/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/ambiguous \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
 300 5
```
-->



luster

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.


```bash
# for DataDir in $(ls -d demulti_dna/*/*/*/*/* | grep -v 'all_loci'); do
  for DataDir in $(ls -d demulti_dna/plate2c/*/*/*/* | grep -v 'multiplex'); do
    # Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    # while [ $Jobs -gt 5 ]; do
    #   sleep 10s
    #   printf "."
    #   Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    # done
    # printf "\n"
    Locus=$(echo $DataDir | cut -f3 -d '/' | sed 's/og/OG/g' | sed 's/TEF1a/TEF/g')
    Prefix=$(echo $DataDir | cut -f3,4,5,6 -d '/' | sed 's&/&_&g' | sed "s/og/OG/g" | sed 's/pathOG/pathog/g' | sed "s/TEF1a/TEF/g")_${Locus}
    # WorkDir=/data2/scratch2/armita/fusarium_ampseq
    R1=$(ls $DataDir/*_${Locus}.fq | head -n1 | tail -n1)
    R2=$(ls $DataDir/*_${Locus}.fq | head -n2 | tail -n1)
    echo $DataDir
    echo $(basename $R1)
    echo $(basename $R2)
    OutDir="processed_dna/"$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev | sed "s/og/OG/g" | sed 's/pathOG/pathog/g' | sed 's/TEF1a/TEF/g')
    # echo $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
  done

for DataDir in $(ls -d demulti_dna/plate2c/*/*/*/* | grep 'multiplex'); do
  for R1 in $(ls $DataDir/*_R1*.fq | grep -v 'ambiguous'); do
    R2=$(echo $R1 | sed 's/_R1/_R2/g')
    echo $(basename $R1)
    echo $(basename $R2)
    Locus=$(echo $R1 | rev | cut -f1 -d '_' | rev | tr -d '.fq')
    Prefix=$(echo $DataDir | cut -f3,4,5,6 -d '/' | sed 's&/&_&g')"_${Locus}"
    OutDir="processed_dna/"$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
    # Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    # while [ $Jobs -gt 1 ]; do
    # sleep 10s
    # printf "."
    # Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    # done
    # printf "\n"
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
  done
done
```


 Summise reads merged:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\tITS\tTef\tSIX13\tT4\tT6\tT8\tAmbiguous\n" > processed_dna/plate2c/merged_summary.tsv
for RunDir in $(ls -d demulti_dna/plate2c/*/*/*/*); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev | sed 's/TEF1a/TEF/g' | sed 's/og/OG/' | sed 's/pathOG/pathog/g')
  Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName="processed_dna/plate2c/$Locus/$Pool/$Dilution/merged/${Locus}_${Pool}_${Dilution}_${Rep}*"
  # RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  ItsLines=$(cat ${RunName}*_*ITS_prefilter.fa | grep '>'| wc -l)
  TefLines=$(cat ${RunName}*_*TEF_prefilter.fa | grep '>'| wc -l)
  Six13=$(ls $RunDir/*_R1_*SIX13.fq)
  Six13Lines=$(cat ${RunName}*_*SIX13_prefilter.fa | grep '>'| wc -l)
  T4Lines=$(cat ${RunName}*_*OG13890_prefilter.fa | grep '>'| wc -l)
  T6Lines=$(cat ${RunName}*_*OG4952_prefilter.fa | grep '>'| wc -l)
  T8Lines=$(cat ${RunName}*_*OG13397_prefilter.fa | grep '>'| wc -l)
  printf "${Locus}_${Pool}_${Dilution}_${Rep}\t$Locus\t$Pool\t$Dilution\t$Rep\t$ItsLines\t$TefLines\t$Six13Lines\t$T4Lines\t$T6Lines\t$T8Lines\n"
done >> processed_dna/plate2c/merged_summary.tsv
ls $PWD/processed_dna/plate2c/merged_summary.tsv
```


Summise reads filtered:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\tITS\tTef\tSIX13\tT4\tT6\tT8\tAmbiguous\n" > processed_dna/plate2c/filtered_summary.tsv
for RunDir in $(ls -d demulti_dna/plate2c/*/*/*/*); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev | sed 's/TEF1a/TEF/g' | sed 's/og/OG/' | sed 's/pathOG/pathog/g')
  Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName="processed_dna/plate2c/$Locus/$Pool/$Dilution/filtered/${Locus}_${Pool}_${Dilution}_${Rep}*"
  ItsLines=$(cat ${RunName}*_*ITS.filtered.fa | grep '>'| wc -l)
  TefLines=$(cat ${RunName}*_*TEF.filtered.fa | grep '>'| wc -l)
  Six13Lines=$(cat ${RunName}*_*SIX13.filtered.fa | grep '>'| wc -l)
  T4Lines=$(cat ${RunName}*_*OG13890.filtered.fa | grep '>'| wc -l)
  T6Lines=$(cat ${RunName}*_*OG4952.filtered.fa | grep '>'| wc -l)
  T8Lines=$(cat ${RunName}*_*OG13397.filtered.fa | grep '>'| wc -l)
  printf "${Locus}_${Pool}_${Dilution}_${Rep}\t$Locus\t$Pool\t$Dilution\t$Rep\t$ItsLines\t$TefLines\t$Six13Lines\t$T4Lines\t$T6Lines\t$T8Lines\n"
done >> processed_dna/plate2c/filtered_summary.tsv
ls $PWD/processed_dna/plate2c/filtered_summary.tsv
```



## OTU assignment
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting

 * All read files associated with the each locus in the project are concatenated
 * Clustering run on all the data for this locus within the project
 * Quantification can then be performed for each treatment against the total set

```bash
  # Concatenate files
  for Locus in ITS TEF SIX13 OG13397 OG13890 OG4952; do
    OutDir=clustering/plate2c/$Locus
    mkdir -p $OutDir
    cat processed_dna/plate2c/$Locus/*/*/filtered/*.filtered.fa > $OutDir/${Locus}_concatenated.temp.fa
    cat processed_dna/plate2c/multiplex_*/*/*/filtered/*_${Locus}.filtered.fa >> $OutDir/${Locus}_concatenated.temp.fa
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


Build databases for loci:

ITS
```bash
qlogin
WorkDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd $WorkDir
OutDir=databases/ITS
mkdir -p $OutDir
cd $OutDir
# The utax database
# https://doi.org/10.15156/BIO/587476
wget -N https://files.plutof.ut.ee/doi/A4/31/A43178A5B024B08CA8B2AB66033B7848BF1A7D0D7A19199363BBC197F7A7F30C.zip
# The UNITE+INSD databases
# https://doi.org/10.15156/BIO/587474
# wget -N https://files.plutof.ut.ee/doi/87/C9/87C97D15437BA13125B61403810C6E46D1319B68A0E6E3BC77BFC57C8D0A67A3.zip
unzip *.zip
mv utax_reference_dataset_01.12.2017.fasta utax_fungi_ITS.fasta
ProgDir=/home/deakig/usr/local/bin
# $ProgDir/usearch -makeudb_sintax UNITE_public_01.12.2017.fasta -output unite_fungi_ITS_sintax.udp
$ProgDir/usearch -makeudb_sintax utax_fungi_ITS.fasta -output utax_fungi_ITS_sintax.udp
cd $WorkDir
```

TEF

```bash
  mkdir -p /home/groups/harrisonlab/project_files/fusarium_ampseq/databases/TEF
```

TEF sequences Provided by AndyT and suplimented with additional
sequences from genbank were aligned and the amplicon extracted.
Non-redundant sequences were kept and samples relabeled to reflect
their original species. Taxa were then exported from geneious and
the followin command was used to add taxonomic information to the
fasta file.

```bash
  # String=';tax=d:Fungi,p:Ascomycota,c:Sordariomycetes,o:Hypocreales,f:Nectriaceae,g:Fusarium,s:'
  # cat TEF_amplicon_nr.fasta | sed "s/>.*/&$String&/" | tr -d '>' | sed "s/^F/>F/g" > TEF_amplicon_nr_db.fasta
  cat TEF_amplicon6_nr_db.fasta | sed 's/_(reversed)//g' > TEF_amplicon_nr_db.fasta
```
The database was then copied bak up to the cluster:
(running from local computer)
```bash
scp TEF_amplicon_nr_db.fasta cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/databases/TEF/.
```

```bash
CurDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd databases/TEF
# cp -s TEF_amplicon6_nr_db.fasta TEF_amplicon_nr_db.fasta
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax TEF_amplicon_nr_db.fasta -output TEF_amplicon_nr_db.udp
cd $CurDir
```


SIX13


```bash
mkdir -p /home/groups/harrisonlab/project_files/fusarium_ampseq/databases/SIX13
```

The database for SIX13 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied bak up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp SIX13_amplicon_nr_db.fasta cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/databases/SIX13/.
```

```bash
CurDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd $CurDir/databases/SIX13
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax SIX13_amplicon_nr_db.fasta -output SIX13_amplicon_nr_db.udp
cd $CurDir
```

T2

```bash
mkdir -p /home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG12981
```

The database for T2 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied back up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp T2_amplicon_nr_db.fasta cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG12981/OG12981_amplicon_nr_db.fasta
```

```bash
CurDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd $CurDir/databases/OG12981
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax OG12981_amplicon_nr_db.fasta -output OG12981_amplicon_nr_db.udp
cd $CurDir
```

T4

```bash
mkdir -p /home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG13890
```

The database for T4 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied bak up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp T4_amplicon_nr_db.fasta cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG13890/OG13890_amplicon_nr_db.fasta
```

```bash
CurDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd $CurDir/databases/OG13890
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax OG13890_amplicon_nr_db.fasta -output OG13890_amplicon_nr_db.udp
cd $CurDir
```


T6

```bash
mkdir -p /home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG4952
```

The database for T6 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied bak up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp T6_amplicon_nr_db.fasta cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG4952/OG4952_amplicon_nr_db.fasta
```

```bash
CurDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd $CurDir/databases/OG4952
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax OG4952_amplicon_nr_db.fasta -output OG4952_amplicon_nr_db.udp
cd $CurDir
```


T8

```bash
mkdir -p /home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG13397
```

The database for T8 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied bak up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp T8_amplicon_nr_db.fasta cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/databases/OG13397/OG13397_amplicon_nr_db.fasta
```

```bash
CurDir=/home/groups/harrisonlab/project_files/fusarium_ampseq
cd $CurDir/databases/OG13397
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax OG13397_amplicon_nr_db.fasta -output OG13397_amplicon_nr_db.udp
cd $CurDir
```

```bash
for Locus in ITS TEF SIX13 OG13397 OG13890 OG4952; do
  for Type in zOTUs; do
  # for Type in OTUs zOTUs; do
  OtuFa=$(ls clustering/$Locus/${Locus}_${Type}.fa)
  RefDb=$(ls databases/$Locus/*.udp)
  Prefix=$(basename ${OtuFa%.fa})
  OutDir=$(dirname $OtuFa)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
  # echo "$OtuFa $RefDb $Prefix $OutDir"
  qsub $ProgDir/submit_taxonomy.sh $OtuFa $RefDb $Prefix $OutDir
  done
done
```

Create OTU tables

```bash
for Locus in ITS TEF; do
  for Pool in soil_pathogens; do
    OutDir=quantified/plate2c/$Locus/$Pool
    QueryReads=$OutDir/${Locus}_reads_appended.fa
    mkdir -p $OutDir
    cat processed_dna/plate2c/$Locus/$Pool/*/filtered/*.fa | cut -f1 -d '.' | sed "s/_${Locus}.*/_${Locus}/g" > $QueryReads
    Threshold=210
    Identity=0.97
    for OtuType in zOTUs; do
      RefDb=$(ls clustering/$Locus/${Locus}_${OtuType}_taxa.fa)
      Prefix="${Locus}_${Pool}_${OtuType}"
      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
      qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
    done
  done
done
for Locus in TEF SIX13 OG13397 OG13890 OG4952; do
  for Pool in Fusarium_spp; do
    OutDir=quantified/plate2c/$Locus/$Pool
    QueryReads=$OutDir/${Locus}_reads_appended.fa
    mkdir -p $OutDir
    cat processed_dna/plate2c/$Locus/$Pool/*/filtered/*.fa | cut -f1 -d '.' | sed "s/_${Locus}.*/_${Locus}/g" > $QueryReads
    Threshold=210
    Identity=1
    for OtuType in zOTUs; do
      RefDb=$(ls clustering/$Locus/${Locus}_${OtuType}_taxa.fa)
      Prefix="${Locus}_${Pool}_${OtuType}"
      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
      qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
    done
  done
done

for Locus in multiplex_PCR; do
  for Pool in Fusarium_spp; do
    OutDir=quantified/plate2c/$Locus/$Pool
    mkdir -p $OutDir
    Threshold=210
    Identity=1
    for OtuType in zOTUs; do
      for RefDb in $(ls databases/*/*_amplicon_nr_db.fasta | grep -e 'TEF' -e 'SIX13' -e 'OG13890' -e 'OG4952' -e 'OG13397'); do
        Locus2=$(basename $RefDb | cut -f1 -d '_')
        QueryReads=$OutDir/${Locus}_$(basename $RefDb | cut -f1 -d '_')_reads_appended.fa
        cat processed_dna/plate2c/$Locus/$Pool/*/filtered/*$(basename $RefDb | cut -f1 -d '_').filtered.fa | cut -f1 -d '.' | sed "s/${Locus2}.*/${Locus2}/g" > $QueryReads
        Prefix="${Locus}_$(basename $RefDb | cut -f1 -d '_')_${Pool}_${OtuType}"
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
      done
    done
  done
done

for Locus in multiplex_barcoding; do
  for Pool in mixed_pool; do
    OutDir=quantified/plate2c/$Locus/$Pool
    mkdir -p $OutDir
    Threshold=210
    Identity=1
    for OtuType in zOTUs; do
      for RefDb in $(ls databases/*/*_amplicon_nr_db.fasta | grep -e 'TEF' -e 'SIX13' -e 'OG13890' -e 'OG4952' -e 'OG13397'); do
        Locus2=$(basename $RefDb | cut -f1 -d '_')
        QueryReads=$OutDir/${Locus}_$(basename $RefDb | cut -f1 -d '_')_reads_appended.fa
        cat processed_dna/plate2c/$Locus/$Pool/*/filtered/*$(basename $RefDb | cut -f1 -d '_').filtered.fa | cut -f1 -d '.' | sed "s/${Locus2}.*/${Locus2}/g" > $QueryReads
        Prefix="${Locus}_$(basename $RefDb | cut -f1 -d '_')_${Pool}_${OtuType}"
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
      done
    done
  done
done
```

```bash
  rm quantified/plate2c/*/*/*.fa
  rm quantified/plate2c/*/*/*_hits.out
```

The output directories were then downloaded to my local computer, where plots
were generated using Rstudio.

on Local machine:
```bash
cd /Users/armita/Downloads/AHDB/ampseq
scp -r cluster:/home/groups/harrisonlab/project_files/fusarium_ampseq/quantified plate2c/.
```

R commands are documented in:
fusarium_ampseq/scripts/plot_OTUs.r

This was run for normalised and unormalised data and using OTU and zOTU data.



```bash
cd AHDB_new/plate2c
for File in $(ls quantified/plate2c/*/*/*_table.txt | grep -v 'multiplex'); do
Primers=$(echo $File | cut -f4 -d '/')
Prefix=$(basename $File | sed 's/_table.txt//g')
OutDir=$(dirname $File)
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
# $ProgDir/plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix --threshold 210
$ProgDir/plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix --threshold 50
done

for File in $(ls quantified/plate2c/*/*/*_table.txt | grep 'multiplex_barcoding'); do
Primers=$(echo $File | cut -f4 -d '/')
Prefix=$(basename $File | sed 's/_table.txt//g')
OutDir=$(dirname $File)
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
# $ProgDir/plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix --threshold 210
$ProgDir/plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix --threshold 50
done

for File in $(ls quantified/plate2c/*/*/*_table.txt | grep 'multiplex_PCR'); do
Primers=$(echo $File | cut -f4 -d '/')
Prefix=$(basename $File | sed 's/_table.txt//g')
OutDir=$(dirname $File)
echo $Prefix
ProgDir=~/cluster_mount/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
# $ProgDir/plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix --threshold 210
$ProgDir/plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix --threshold 50
done

```
