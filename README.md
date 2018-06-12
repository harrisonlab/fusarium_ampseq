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

RUN=.
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
```



```