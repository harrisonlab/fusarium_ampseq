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

# Analysis of individual loci

## Demultiplexing

This script demultiplexs mixed (e.g. ITS and 16S) libraries based on the primer sequence. Any sequence which has mismatches is written to ambiguous.fq (f & r seperately). Primer sequences
are detailed in the submission wrapper.
*Note* Regex are used to describe degenerate bases in the primer.

Run below to demultiplex:

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
    qsub $ProgDir/submit_demulti.sh $R1 $R2 ${OutDir}
  done
```
 Summise reads demultiplexed:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\tITS\tTef\tSIX13\tT2\tT4\tAmbiguous\n" > demulti_dna/demultiplex_summary.tsv
for RunDir in $(ls -d demulti_dna/*/*/*/*); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  ITS=$(ls $RunDir/*ITS.fq)
  ItsLines=$(cat $ITS | wc -l)
  TEF=$(ls $RunDir/*TEF.fq)
  TefLines=$(cat $TEF | wc -l)
  Six13=$(ls $RunDir/*SIX13.fq)
  Six13Lines=$(cat $Six13 | wc -l)
  T2=$(ls $RunDir/*OG12981.fq)
  T2Lines=$(cat $T2 | wc -l)
  T4=$(ls $RunDir/*OG13890.fq)
  T4Lines=$(cat $T4 | wc -l)
  Ambiguous=$(ls $RunDir/*ambiguous.fq)
  AmbLines=$(cat $Ambiguous | wc -l)
  printf "$RunName\t$Locus\t$Pool\t$Dilution\t$Rep\t$ItsLines\t$TefLines\t$Six13Lines\t$T2Lines\t$T4Lines\t$AmbLines\n"
done >> demulti_dna/demultiplex_summary.tsv
```

From this data thresholding values of cross-contamination as a result of illumina
adapter read hopping were determined. For each row in the dataset, the most
abundant locus was assumed to be the target locus and reads attributed to other
loci assumed to be contaminant reads. Reads from another experiment with the same
locus could be contaminating the sample. The 2nd most abundant locus was
identified for each run (the highest contaminant locus) and the maximum value
identified accross the entire plate. A threshold for a minimum abundance of reads
attributed to an OTU was set at this value.

In the case of our plate, this was:
```bash
Threshold=206
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


## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.


```bash
for DataDir in $(ls -d demulti_dna/*/*/*/* | grep -v 'all_loci'); do
  # Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
  # while [ $Jobs -gt 1 ]; do
  # sleep 1m
  # printf "."
  # Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
  # done
  printf "\n"
  Locus=$(echo $DataDir | cut -f2 -d '/')
  Prefix=$(echo $DataDir | cut -f2,3,4,5 -d '/' | sed 's&/&_&g')
  WorkDir=/data/scratch/armita/fusarium_ampseq
  R1=$(ls $DataDir/*_${Locus}.fq | head -n1 | tail -n1)
  R2=$(ls $DataDir/*_${Locus}.fq | head -n2 | tail -n1)
  echo $DataDir
  echo $(basename $R1)
  echo $(basename $R2)
  OutDir="processed_dna/"$(echo $R1 | rev | cut -f2,3,4,5 -d '/' | rev)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
  qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
done
```


## OTU assignment
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting

 * All read files associated with the each locus in the project are concatenated
 * Clustering run on all the data for this locus within the project
 * Quantification can then be performed for each treatment against the total set

```bash
  # Concatenate files
  # for Locus in "ITS" "TEF"; do
  for Locus in SIX13 OG12981 OG13890; do
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

TEF

TEF sequences Provided by AndyT and suplimented with additional
sequences from genbank were aligned and the amplicon extracted.
Non-redundant sequences were kept and samples relabeled to reflect
their original species. Taxa were then exported from geneious and
the followin command was used to add taxonomic information to the
fasta file.

```bash
  String=';tax=d:Fungi,p:Ascomycota,c:Sordariomycetes,o:Hypocreales,f:Nectriaceae,g:Fusarium,s:'
  cat TEF_amplicon_nr.fasta | sed "s/>.*/&$String&/" | tr -d '>' | sed "s/^F/>F/g" > TEF_amplicon_nr_db.fasta
```
The database was then copied bak up to the cluster:
(running from local computer)
```bash
scp TEF_amplicon_nr_db.fasta cluster:/data/scratch/armita/fusarium_ampseq/databases/TEF/.
```

```bash
cd /data/scratch/armita/fusarium_ampseq/databases/TEF
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax TEF_amplicon_nr_db.fasta -output TEF_amplicon_nr_db.udp
cd /data/scratch/armita/fusarium_ampseq
```

SIX13

The database for SIX13 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied bak up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp SIX13_amplicon_nr_db.fasta cluster:/data/scratch/armita/fusarium_ampseq/databases/SIX13/.
```

```bash
cd /data/scratch/armita/fusarium_ampseq/databases/SIX13
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax SIX13_amplicon_nr_db.fasta -output SIX13_amplicon_nr_db.udp
cd /data/scratch/armita/fusarium_ampseq
```


T2

The database for T2 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied bak up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp T2_amplicon_nr_db.fasta cluster:/data/scratch/armita/fusarium_ampseq/databases/OG12981/OG12981_amplicon_nr_db.fasta
```

```bash
cd /data/scratch/armita/fusarium_ampseq/databases/OG12981
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax OG12981_amplicon_nr_db.fasta -output OG12981_amplicon_nr_db.udp
cd /data/scratch/armita/fusarium_ampseq
```

T4

The database for T4 was made using the same approach
as TEF. Sequences originated from reference genomes.
The database was then copied bak up to the cluster:
(running from local computer)

The database was copied to the cluster:
(running from local computer)
```bash
scp T4_amplicon_nr_db.fasta cluster:/data/scratch/armita/fusarium_ampseq/databases/OG13890/OG13890_amplicon_nr_db.fasta
```

```bash
cd /data/scratch/armita/fusarium_ampseq/databases/OG13890
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax OG13890_amplicon_nr_db.fasta -output OG13890_amplicon_nr_db.udp
cd /data/scratch/armita/fusarium_ampseq
```

```bash
for Locus in OG12981; do
# for Locus in ITS TEF SIX13 OG12981 OG13890; do
  for Type in OTUs zOTUs; do
  OtuFa=$(ls clustering/$Locus/${Locus}_${Type}.fa)
  RefDb=$(ls databases/$Locus/*.udp)
  Prefix=$(basename ${OtuFa%.fa})
  OutDir=$(dirname $OtuFa)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
  qsub $ProgDir/submit_taxonomy.sh $OtuFa $RefDb $Prefix $OutDir
  done
done
```

Create OTU tables

```bash
# for Locus in ITS TEF SIX13 OG12981 OG13890; do
for Locus in OG12981; do
  for Pool in soil_pathogens Fusarium_spp; do
  # for Pool in Fusarium_spp; do
    OutDir=quantified/$Locus/$Pool
    mkdir -p $OutDir
    cat processed_dna/$Locus/$Pool/*/*/merged/*.fa | cut -f1 -d '.' > $OutDir/${Locus}_reads_appended.fa
    QueryReads=$(ls $OutDir/${Locus}_reads_appended.fa)
    # for OtuType in OTUs zOTUs; do
    for OtuType in zOTUs; do
      RefDb=$(ls clustering/$Locus/${Locus}_${OtuType}_taxa.fa)
      Prefix=$Locus
      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
      qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir
    done
  done
done
```

```bash
# for Locus in ITS TEF SIX13 OG12981 OG13890; do
for Locus in OG12981; do
  for Pool in soil_pathogens Fusarium_spp; do
  # for Pool in soil_pathogens; do
    for OutDir in $(ls -d quantified/$Locus/$Pool); do
      # for OtuType in OTUs zOTUs; do
      for OtuType in zOTUs; do
        Prefix=$Locus
        # Plot by species
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        $ProgDir/plot_OTUs.r --OTU_table $OutDir/${Prefix}_${OtuType}_table_by_spp.txt --prefix $OutDir/${Prefix}_${OtuType}_table_by_spp
        $ProgDir/plot_OTUs.r --OTU_table $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded.txt --prefix $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded
        $ProgDir/plot_OTUs.r --OTU_table $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded_norm.txt --prefix $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded_norm
      done
    done
  done
done
```

```bash
rm quantified/*/*/*_reads_appended.fa
rm quantified/*/*/*_hits.out
```



The output directories were then downloaded to my local computer, where plots
were generated using Rstudio.

on Local machine:
```bash
cd /Users/armita/Downloads/AHDB/ampseq
scp -r cluster:/data/scratch/armita/fusarium_ampseq/quantified/ITS/soil_pathogens .
```

R commands are documented in:
fusarium_ampseq/scripts/plot_OTUs.r

This was run for normalised and unormalised data and using OTU and zOTU data.



# Analysis of pooled loci

## Demultiplexing

This script demultiplexs mixed (e.g. ITS and 16S) libraries based on the primer sequence. Any sequence which has mismatches is written to ambiguous.fq (f & r seperately). Primer sequences
are detailed in the submission wrapper.
*Note* Regex are used to describe degenerate bases in the primer.

Run below to demultiplex:

```bash
  for DataDir in $(ls -d raw_dna/paired/*/*/*/* | grep 'all_loci'); do
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
    qsub $ProgDir/submit_demulti.sh $R1 $R2 ${OutDir}
  done
```
 Summise reads demultiplexed:
```bash
printf "RunName\tLocus\tPool\tDilution\tRep\tITS\tTef\tSIX13\tT2\tT4\tAmbiguous\n" > demulti_dna/demultiplex_summary.tsv
for RunDir in $(ls -d demulti_dna/*/*/*/* | grep 'all_loci'); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Pool=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  ITS=$(ls $RunDir/*ITS.fq)
  ItsLines=$(cat $ITS | wc -l)
  TEF=$(ls $RunDir/*TEF.fq)
  TefLines=$(cat $TEF | wc -l)
  Six13=$(ls $RunDir/*SIX13.fq)
  Six13Lines=$(cat $Six13 | wc -l)
  T2=$(ls $RunDir/*OG12981.fq)
  T2Lines=$(cat $T2 | wc -l)
  T4=$(ls $RunDir/*OG13890.fq)
  T4Lines=$(cat $T4 | wc -l)
  Ambiguous=$(ls $RunDir/*ambiguous.fq)
  AmbLines=$(cat $Ambiguous | wc -l)
  printf "$RunName\t$Locus\t$Pool\t$Dilution\t$Rep\t$ItsLines\t$TefLines\t$Six13Lines\t$T2Lines\t$T4Lines\t$AmbLines\n"
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
```
-->


## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.


```bash
for DataDir in $(ls -d demulti_dna/*/*/*/* | grep 'all_loci'); do
  for R1 in $(ls $DataDir/*_R1_*.fq | grep -v 'ambiguous'); do
    R2=$(echo $R1 | sed 's/_R1_/_R2_/g')
    Locus=$(echo ${R1%.fq} | rev | cut -f1 -d '_' | rev)
    Prefix=$(echo $DataDir | cut -f2,3,4,5 -d '/' | sed 's&/&_&g')_${Locus}
    # Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    # while [ $Jobs -gt 1 ]; do
    # sleep 1m
    # printf "."
    # Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    # done
    # printf "\n"
    WorkDir=/data/scratch/armita/fusarium_ampseq
    echo $DataDir
    echo $(basename $R1)
    echo $(basename $R2)
    OutDir="processed_dna/"$(echo $R1 | rev | cut -f2,3,4,5 -d '/' | rev)"/${Locus}"
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
  done
done
```


## OTU assignment
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. Greg has written his own scripts to do the dereplication and sorting

 * All read files associated with the each locus in the project are concatenated
 * Clustering run on all the data for this locus within the project
 * Quantification can then be performed for each treatment against the total set

 * This may not be the best approach as zOTUs are identified by separating out real and error reads. Therefore this will be sensitive to taxa abundance within the dataset. For example, we would not want a taxon that appears in only one dataset to be identified as error-reads by abundance of a similar taxon in all other datasets.

```bash
  # Concatenate files
  for Locus in ITS TEF SIX13 OG12981 OG13890; do
    OutDir=clustering/all_loci/$Locus
    mkdir -p $OutDir
    cat processed_dna/all_loci/*/*/*/$Locus/filtered/*.filtered.fa > $OutDir/${Locus}_concatenated.temp.fa
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
* See commands above


```bash
for Locus in ITS TEF SIX13 OG12981 OG13890; do
  for Type in OTUs zOTUs; do
  OtuFa=$(ls clustering/all_loci/$Locus/${Locus}_${Type}.fa)
  RefDb=$(ls databases/$Locus/*.udp)
  Prefix=$(basename ${OtuFa%.fa})
  OutDir=$(dirname $OtuFa)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
  qsub $ProgDir/submit_taxonomy.sh $OtuFa $RefDb $Prefix $OutDir
  done
done
```

Create OTU tables

```bash
for Locus in ITS TEF SIX13 OG12981 OG13890; do
  for Pool in mixed_pool Fusarium_spp; do
    OutDir=quantified/all_loci/$Pool/$Locus
    mkdir -p $OutDir
    cat processed_dna/all_loci/$Pool/*/*/$Locus/merged/*.fa | cut -f1 -d '.' > $OutDir/${Locus}_reads_appended.fa
    QueryReads=$(ls $OutDir/${Locus}_reads_appended.fa)
    for OtuType in OTUs zOTUs; do
      RefDb=$(ls clustering/$Locus/${Locus}_${OtuType}_taxa.fa)
      Prefix=$Locus
      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
      qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir
    done
  done
done
```


```bash
for Locus in ITS TEF SIX13 OG12981 OG13890; do
# for Locus in OG12981; do
  for Pool in mixed_pool Fusarium_spp; do
  # for Pool in soil_pathogens; do
    for OutDir in $(ls -d quantified/all_loci/$Pool/$Locus); do
      # for OtuType in OTUs zOTUs; do
      for OtuType in zOTUs; do
        Prefix=$Locus
        # Plot by species
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        $ProgDir/plot_OTUs.r --OTU_table $OutDir/${Prefix}_${OtuType}_table_by_spp.txt --prefix $OutDir/${Prefix}_${OtuType}_table_by_spp
        $ProgDir/plot_OTUs.r --OTU_table $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded.txt --prefix $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded
        $ProgDir/plot_OTUs.r --OTU_table $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded_norm.txt --prefix $OutDir/${Prefix}_${OtuType}_table_by_spp_thresholded_norm
      done
    done
  done
done
```



```bash
rm quantified/*/*/*/*_reads_appended.fa
rm quantified/*/*/*/*_hits.out
```

```bash
for File in $(ls all_loci/*/*/*_table.txt); do
# for File in $(ls all_loci/mixed_pool/*/*_table.txt); do
Locus="all_loci"
Pool=$(echo $File | cut -f3 -d '/')
OtuType=$(basename $File | cut -f2 -d '_')
Primers=$(echo $File | cut -f4 -d '/')
Prefix=${Locus}_${Pool}_${OtuType}_${Primers}
OutDir=$(dirname $File)
# ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
./plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix --threshold 206
done
for File in $(ls all_loci/*/*/*_table_norm.txt); do
# for File in $(ls all_loci/mixed_pool/*/*_table.txt); do
Locus="all_loci"
Pool=$(echo $File | cut -f3 -d '/')
OtuType=$(basename $File | cut -f2 -d '_')
Primers=$(echo $File | cut -f4 -d '/')
Prefix="${Locus}_${Pool}_${OtuType}_${Primers}_norm"
OutDir=$(dirname $File)
# ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
./plot_OTUs_local.r --OTU_table $File --prefix $OutDir/$Prefix
done
```

The output directories were then downloaded to my local computer, where plots
were generated using Rstudio.

on Local machine:
```bash
cd /Users/armita/Downloads/AHDB/ampseq
scp -r cluster:/data/scratch/armita/fusarium_ampseq/quantified/ITS/Fusarium_spp/*.txt .
scp -r cluster:/data/scratch/armita/fusarium_ampseq/quantified/ITS/mixed_pool/*.txt .
```

R commands are documented in:
fusarium_ampseq/scripts/plot_OTUs.r

This was run for normalised and unormalised data and using OTU and zOTU data.
