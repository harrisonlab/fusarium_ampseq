# fusarium_ampseq
Commands used in analysis of amplicon sequence data of Fusarium infested soil samples and artificial mixes.


# Data transfer

Samples sequenced were:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180427_M04465_0077_000000000-BM7W8/Data/
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


Raw sequencing data was symbolicly linked to the working directory:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180427_M04465_0077_000000000-BM7W8/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  for ReadF in $(ls $RawDatDir/*_L001_R1_001.fastq.gz | grep -v 'Undetermined'); do
    ReadR=$(echo $ReadF | sed 's/_L001_R1_001.fastq.gz/_L001_R2_001.fastq.gz/g')
    Name=$(basename $ReadF _L001_R1_001.fastq.gz)
    Locus=$(echo $Name | cut -f1 -d '-')
    Sample=$(echo $Name | cut -f2 -d '_' | tr -d 'S')
    Dilution=$(echo $Name | cut -f2 -d '-' | sed "s/.//")
    TechRep=$(echo $Name | cut -f3 -d '-' | cut -f1 -d '_')
    # echo "$Name - $Sample - $Locus - $Dilution - $TechRep"
    # Set species pools
    if [ $Sample -le 24 ]; then
      Pool="soil_pathogens"
    elif [ $Sample -gt 29 ]; then
      Pool="Fusarium_spp"
    else
      Pool="mixed_pool"
    fi
    # Set dilution groups
    if [ $Dilution == 'a' ]; then
      Dilution="equimolar"
    elif [ $Dilution == 'b' ]; then
      Dilution="dilution_10x"
    elif [ $Dilution == 'c' ]; then
      Dilution="dilution_100x"
    elif [ $Dilution == 'P' ] ; then
      Dilution="equimolar"
    fi
    # Set Locus
    if [ $Locus == 'PCR1' ] || [ $Locus == 'PCR2' ] || [ $Locus == 'PPCR1' ]; then
      Locus="all_loci"
    elif [ $Locus == 'T1' ]; then
      Locus="SIX13"
    elif [ $Locus == 'T2' ]; then
      Locus="OG12981"
      # Locus="OG12981_FoC_Fus2_g17143"
    elif [ $Locus == 'T4' ] ; then
      Locus="OG13890"
      # Locus="OG13890_FoC_Fus2_PGN_07490"
    fi
    # echo "$Name - $Pool - $Locus - $Dilution - $TechRep"
    # Set technical rep:
    if [ $TechRep == 'Pool' ]; then
      TechRep="pooled_reps"
    fi
    echo "$Locus - $Pool - $Dilution - $TechRep"
    OutDir=/data/scratch/armita/fusarium_ampseq/raw_dna/paired/$Locus/$Pool/$Dilution/$TechRep
    mkdir -p $OutDir/F
    cp -s $ReadF $OutDir/F/.
    mkdir -p $OutDir/R
    cp -s $ReadR $OutDir/R/.
  done
```


#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
  cd /data/scratch/armita/fusarium_ampseq
  for RawData in $(ls raw_dna/paired/*/*/*/*/*/*.fastq.gz); do
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
for StrainPath in $(ls -d raw_dna/paired/*/*/*/*); do
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
for F_read in $(ls qc_dna/paired/*/*/*/*/F/*fq.gz); do
  Locus=$(echo $F_read | cut -f3 -d '/')
  Pool=$(echo $F_read | cut -f4 -d '/')
  Dilution=$(echo $F_read | cut -f5 -d '/')
  TechRep=$(echo $F_read | cut -f6 -d '/')
  ReadCount=$(cat $F_read | gunzip -cf | awk '{s++}END{print s/4}')
  printf "$Locus\t$Pool\t$Dilution\t$TechRep\t$ReadCount\n"
done > qc_dna/reads_per_sample.tsv
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

# Common workflow

## Set up project folders

If the project has multiple sequencing runs, RUN should be set to location where files are to be stored.

```shell
PROJECT_FOLDER=/home/groups/harrisonlab/project_files/fusarium_ampseq
mkdir -p $PROJECT_FOLDER
ln -s $MBPL $PROJECT_FOLDER/metabarcoding_pipeline

RUN=test
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
mkdir $PROJECT_FOLDER/data/$RUN/quality
mkdir $PROJECT_FOLDER/data/$RUN/ambiguous

for s in BAC FUN OO NEM; do
  mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
  mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
  mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done
```
<!--
Copy raw Fastq files (or link to) to $PROJECT_FOLDER/data/$RUN/fastq

## Decompress files (not required, unless using none multiplexed data)

Append 2 to decompress sym links

```shell
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*.gz; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c unzip $FILE 2
done
```

## QC
Qualtiy checking with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```shell
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
done
``` -->




## Demultiplexing

Script demulti.pl demultiplexs mixed (e.g. ITS and 16S) libraries based on the primer sequence. Number of acceptable mismatches in the primer sequence can be specified (0 by default). Any sequence which has too many mismatches, or none mathching primers is written to ambiguous.fq (f & r seperately). The script accepts multiple primer pairs.

<table>
Possible primers:
<tr><td><td>Forward<td>Reverse</tr>
<tr><td>Bacteria<td>CCTACGGGNGGCWGCAG<td>GACTACHVGGGTATCTAATCC</tr>
<tr><td>Fungi<td>CTTGGTCATTTAGAGGAAGTAA<td>ATATGCTTAAGTTCAGCGGG</tr>
<tr><td>Oomycete<td>GAAGGTGAAGTCGTAACAAGG<td>AGCGTTCTTCATCGATGTGC</tr>
<tr><td>Nematode<td>CGCGAATRGCTCATTACAACAGC<td>GGCGGTATCTGATCGCC</tr>
</table>

If the primers are unknown, running something like the below should give a good indication of what they are (run from folder of decompressed fastq files). It will also give a good indication of how many mismatches (if any) to use for demulti.pl.
```shell
sed -n '2~4p' $(ls|head -n1)|grep -x "[ATCG]\+"|cut -c-16|sort|uniq| \
tee zzexpressions.txt|xargs -I%  grep -c "^%" $(ls|head -n1) >zzcounts.txt
```

CCGGCCTACTGGTTTCTGAGCGCAGCACAAGTCGCGCCGGCCTACTGGTTTCTGAGCGCAGCACAAGTCGCGCCGGCCTACTGGTTTCTGAGCGCAGCACAAGTCGCG
```
ITS primer with illumina adapter F
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGAATCATCGAATCTTTGAACGC
ITS primer with illumina adapter R
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCGCTTATTGATATGCTTAARTTCAG
```

Run below to demultiplex

Note - this step may be problematic for future runs as we have
  * multiple F or R primers for some loci
  * degenerate primer site for some primers
As such, I may need to write a new demultiplexing script.

```bash
  for DataDir in $(ls -d raw_dna/paired/*/*/*/*); do
    Jobs=$(qstat | grep 'submit_dem' | grep 'qw'| wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 1m
    printf "."
    Jobs=$(qstat | grep 'submit_dem' | grep 'qw'| wc -l)
    done
    printf "\n"
    WorkDir=/data/scratch/armita/fusarium_ampseq
    R1=$(ls $DataDir/F/*.fastq.gz)
    R2=$(ls $DataDir/R/*.fastq.gz)
    echo $DataDir
    echo $(basename $R1)
    echo $(basename $R2)
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    OutDir=demulti_dna/$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
    # echo $OutDir
    qsub $ProgDir/submit_demulti.sh $R1 $R2 ${OutDir}
  done
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
``` -->

# Oomycete ITS workflow


## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size min_join_overlap max_errrors

```shell

# PROJECT_FOLDER=/home/groups/harrisonlab/project_files/fusarium_ampseq
# ProgDir=/home/armita/git_repos/emr_repos/scripts/Metabarcoding_pipeline/scripts
# $ProgDir/PIPELINE.sh -c OOpre \
#   "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
#   $PROJECT_FOLDER/data/$RUN/$SSU \
#   $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
#   $MINL $MINOVER $QUAL $FPL $RPL

# R1=$(ls demulti_dna/ITS/soil_pathogens/equimolar/1/ITS-1a-1_S1_L001_R1_001_ITS.fq)
# R2=$(ls demulti_dna/ITS/soil_pathogens/equimolar/1/ITS-1a-1_S1_L001_R2_001_ITS.fq)
for DataDir in $(ls -d demulti_dna/*/*/*/* | grep 'ITS'); do
  Jobs=$(qstat | grep 'sub_process' | grep 'qw'| wc -l)
  while [ $Jobs -gt 1 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_process' | grep 'qw'| wc -l)
  done
  printf "\n"
  Locus=$(echo $DataDir | cut -f2 -d '/')
  Prefix=$(echo $DataDir | cut -f2,3,4,5 -d '/' | sed 's&/&_&g')
  WorkDir=/data/scratch/armita/fusarium_ampseq
  R1=$(ls $DataDir/*_${Locus}.fq | head -n1 | tail -n1)
  R2=$(ls $DataDir/*_${Locus}.fq | head -n2 | tail -n1)
  echo $DataDir
  echo $(basename $R1)
  echo $(basename $R2)
# R1=$(ls demulti_dna/TEF/Fusarium_spp/equimolar/1/TEF-2a-1_S65_L001_R1_001_TEF.fq)
# R2=$(ls demulti_dna/TEF/Fusarium_spp/equimolar/1/TEF-2a-1_S65_L001_R2_001_TEF.fq)
  OutDir="processed_dna/"$(echo $R1 | rev | cut -f2,3,4,5 -d '/' | rev)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
  qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
done

#
# #S=$(echo $f|awk -F"_" -v D=$(echo $LOC|awk -F"/" '{print($(NF-3))}') '{print $2"D"D}')
# S=$(echo $f|awk -F"/" '{print $NF'}|awk -F"_" {'print $1,$2'} OFS="_")
#
# OUTFILE=$1
# OUTDIR=$1
# ADAPTERS=$1
# MINL=$1
# MAXDIFF=$1
# QUAL=$1
# FPL=$1
# RPL=$1
# SCRIPT_DIR=$1
# qsub $SCRIPT_DIR/submit_16Spre_v2.sh $R1 $R2 $OUTFILE $S $@ $SCRIPT_DIR
```


## OTU assignment
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting

 * All read files associated with the each locus in the project are concatenated
 * Clustering run on all the data for this locus within the project
 * Quantification can then be performed for each treatment against the total set

```bash
  # Concatenate files
  for Locus in "ITS"; do
    OutDir=clustering/$Locus
    mkdir -p $OutDir
    cat processed_dna/$Locus/*/*/*/filtered/*.filtered.fa > $OutDir/${Locus}_concatenated.temp.fa
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
WorkDir=/data/scratch/armita/fusarium_ampseq
cd $WorkDir
OutDir=databases/ITS
mkdir -p $OutDir
cd $OutDir
# https://doi.org/10.15156/BIO/587476
wget -N https://files.plutof.ut.ee/doi/A4/31/A43178A5B024B08CA8B2AB66033B7848BF1A7D0D7A19199363BBC197F7A7F30C.zip
unzip *.zip
mv utax_reference_dataset_01.12.2017.fasta utax_fungi_ITS.fasta
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax utax_fungi_ITS.fasta -output utax_fungi_ITS_sintax.udp
cd $WorkDir
```

```bash
# $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU sintax
# qsub $SCRIPT_DIR/submit_taxonomy.sh $SCRIPT_DIR $@
OTUs=clustering/ITS/ITS_OTUs.fa
RefDb=$(ls databases/ITS/utax_fungi_ITS_sintax.udp)
Prefix=$(echo $OTUs | cut -f2 -d '|')"_OTUs"
OutDir=$(dirname $OTUs)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
qsub $ProgDir/submit_taxonomy.sh $OTUs $RefDB $Prefix $OutDir

```

Create OTU tables
```bash
# $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN $SSU 21 20
```
```bash
Primers=primers.fa
printf \
">ITS_F
GTGAATCATCGAATCTTTGAACGC
>ITS_R
CCGCTTATTGATATGCTTAARTTCAG
>TEF_F
GGTCACTTGATCTACCAGTGCG
>TEF_R
CCCARGCGTACTTGAAGRAAC
>SIX_F
GCTACTCAAAGTCGTGGACGAG
>SIX_R
GGCAATATATTCCGTCCATTCTTGG
>FOCg17143_F
CACTTCCTCACTTACTTTACCACTCC
>FOCg17143_F
GTCATCGCAATCGCCKTCCG
>orthogroup_13890_F
GCTGTCTTATCACTTATCAGCCTTG
>orthogroup_13890_R
CGGTCTGATTTGGTGTCCAGTCG" \
> $Primers

  FilteredReads=$(ls processed_dna/ITS/soil_pathogens/equimolar/1/filtered/ITS_soil_pathogens_equimolar_1.filtered.fa)
  UnfilteredReads=$(ls processed_dna/ITS/soil_pathogens/equimolar/1/unfiltered/ITS_soil_pathogens_equimolar_1.unfiltered.fastq)
  OutDir=$(basename $FilteredReads)
	Prefix=$(echo $DataDir | cut -f2,3,4,5 -d '/' | sed 's&/&_&g')
  Locus=$(echo $DataDir | cut -f2 -d '/')
  FPL=$(cat $Primers | grep -A1 "${Locus}_F" | tail -n1 | wc -c)
  RPL=$(cat $Primers | grep -A1 "${Locus}_R" | tail -n1 | wc -c)

	EP=both

	JOBNAME=OTU_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

	cd $OutDir
	dir=`mktemp -d -p $OutDir`
	find $UNFILTDIR -name '*.fastq' >$dir/files.txt
	TASKS=$(wc -l $dir/files.txt|awk -F" " '{print $1}')
	qsub -N ${JOBNAME}_1 -t 1-$TASKS:1 $SCRIPT_DIR/submit_fastq_fasta.sh $dir/files.txt $dir $FPL $RPL $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_1 -N ${JOBNAME}_2 $SCRIPT_DIR/submit_cat_files.sh $dir $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_2 -N ${JOBNAME}_3 $SCRIPT_DIR/submit_global_search.sh $dir/t1.fa $OutDir $PREFIX otus plus
	qsub -hold_jid ${JOBNAME}_2 -N ${JOBNAME}_4 $SCRIPT_DIR/submit_global_search.sh $dir/t1.fa $OutDir $PREFIX zotus plus

	if $R2; then
		find $UNFILTDIR -name '*.r2.*' >$dir/R2.files.txt			
		TASKS=$(wc -l $dir/R2.files.txt|awk -F" " '{print $1}')
		qsub -hold_jid ${JOBNAME}_3 -N ${JOBNAME}_4 -t 1-$TASKS:1 $SCRIPT_DIR/submit_search_hits.sh $dir/files.txt $dir $OutDir/$PREFIX.hits.out $SCRIPT_DIR
		qsub -hold_jid ${JOBNAME}_4 -N ${JOBNAME}_5 $SCRIPT_DIR/submit_global_search.sh $dir/t3.fa $OutDir $PREFIX otus both
		qsub -hold_jid ${JOBNAME}_4 -N ${JOBNAME}_6 $SCRIPT_DIR/submit_global_search.sh $dir/t3.fa $OutDir $PREFIX zotus both
		qsub -hold_jid ${JOBNAME}_5,${JOBNAME}_6 $SCRIPT_DIR/submit_tidy.sh $dir $PREFIX.hits.out ${PREFIX}2.hits.out OTU_*_1.* OTU_*_4.*
	else
		qsub -hold_jid ${JOBNAME}_3,${JOBNAME}_4 $SCRIPT_DIR/submit_tidy.sh $dir $PREFIX.hits.out OTU_*_1.*
	fi  

	exit 1
	;;
```
