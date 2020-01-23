# fusarium_ampseq

Commands used in analysis of amplicon sequence data of Fusarium infested soil samples and artificial mixes as part of AHDB project 2019-2020.

All of this work was done in the directory:
```bash
  mkdir -p /projects/fusarium_ampseq
  cd /projects/fusarium_ampseq
```
<!--
# Identify presence of SIX5 in Fusarium genomes

The locus SIX5 was selected for identification. A qPCR diagnostic has been
developed for locus, which has been shown to be specific, but a new reverse
primer has been developed. As such, the presence of SIX5 was searched accross
sequenced Fusarium genomes.

```bash


``` -->


# Data transfer

Samples sequenced were:

```bash
  RawDat=$(ls -d /archives/2019_niabemr_miseq/RAW/191220_M04465_0110_000000000-CPB7N/Data/Intensities/BaseCalls /archives/2019_niabemr_miseq/RAW/200103_M04465_0111_000000000-CVMTT/Data/Intensities/BaseCalls)
  ls $RawDatDir
```

Sample names referred to the following:

Loci:
16S, ITS, TEF, SIX5 (SIX13), T4 (OG13890), T6 (OG4952) with combined 16S and ITS samples run using a single barcode and TEF using a single barcode with multiplexed SIX5, T4, T6. Four Soils were run accross the two plates with pathogen inoculations at 4 different concentrations for each soil. Additionally, soils from previous analyses were included on both plates.


## Sample assignment

```bash
WorkDir=/projects/fusarium_ampseq

for FileF in $(ls /archives/2019_niabemr_miseq/RAW/191220_M04465_0110_000000000-CPB7N/Data/Intensities/BaseCalls/*_L001_R1_001.fastq.gz /archives/2019_niabemr_miseq/RAW/200103_M04465_0111_000000000-CVMTT/Data/Intensities/BaseCalls/*_L001_R1_001.fastq.gz | grep -v 'Undetermined'); do
  FileR=$(echo $FileF | sed "s/_R1_/_R2_/g")
  Date=$(echo $FileF | cut -f5 -d '/'| cut -f1 -d '_')
  if [[ "$Date" == "191220" ]]; then
    Plate="1"
  else
    Plate="2"
  fi
  Name=$(basename $FileF)
  Field=$(echo $Name | cut -f1 -d '-')
  if [[ "FOC" == "$Field" || "FON" == "$Field" || "FOM" == "$Field" ]]; then
    # echo $Field
    Dilution="NA"
    TechRep=$(echo $Name | cut -f2 -d '-')
    Locus=$(echo $Name | cut -f3- -d '-' | sed "s/_S.*//g")
  else
    Dilution=$(echo $Name | cut -f2 -d '-')
    TechRep=$(echo $Name | cut -f3 -d '-')
    Locus=$(echo $Name | cut -f4- -d '-' | sed "s/_S.*//g")
  fi
  Mix=$(echo "$Locus" | sed "s/16s-ITS/ITS_mix/g" | sed "s/TEF/TEF_mix/g")
  printf "$Name\t$Plate\t$Field\t$Dilution\t$TechRep\t$Locus\n"
  OutDir=$WorkDir/raw_dna/$Plate/paired/$Mix/$Field/$Dilution/$TechRep
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cp -s $FileF $OutDir/F/.
  cp -s $FileR $OutDir/R/.
done > $WorkDir/sample_table.tsv
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
for DataDir in $(ls -d raw_dna/*/paired/*/*/*/*); do
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
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/2020_analysis
OutDir=$(echo $DataDir | sed "s/raw_dna/demulti_dna/g" | sed 's&paired/&&g')
echo $OutDir
qsub $ProgDir/submit_demulti_2020.sh $R1 $R2 ${OutDir}
done
```


 Summise reads demultiplexed:
```bash
printf "RunName\tplate\tLocus\tSoil\tDilution\tRep\t16S\tITS\tTef\tSIX5\tT4\tT6\tAmbiguous\n" > demulti_dna/demultiplex_summary.tsv
for RunDir in $(ls -d demulti_dna/*/*/*/*/*); do
  Plate=$(echo $RunDir | rev | cut -f5 -d '/' | rev)
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Soil=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  T0=$(ls $RunDir/*_R1_*16S.fq)
  T0Lines=$(cat $T0 | awk '{s++}END{print s/4}')
  ITS=$(ls $RunDir/*_R1_*ITS.fq)
  ItsLines=$(cat $ITS | awk '{s++}END{print s/4}')
  TEF=$(ls $RunDir/*_R1_*TEF.fq)
  TefLines=$(cat $TEF | awk '{s++}END{print s/4}')
  Six5=$(ls $RunDir/*_R1_*SIX5.fq)
  Six5Lines=$(cat $Six5 | awk '{s++}END{print s/4}')
  T4=$(ls $RunDir/*_R1_*OG13890.fq)
  T4Lines=$(cat $T4 | awk '{s++}END{print s/4}')
  T6=$(ls $RunDir/*_R1_*OG4952.fq)
  T6Lines=$(cat $T6 | awk '{s++}END{print s/4}')
  Ambiguous=$(ls $RunDir/*_R1_*ambiguous.fq)
  AmbLines=$(cat $Ambiguous | awk '{s++}END{print s/4}')
  printf "$RunName\t$Plate\t$Locus\t$Soil\t$Dilution\t$Rep\t$T0Lines\t$ItsLines\t$TefLines\t$Six5Lines\t$T4Lines\t$T6Lines\t$AmbLines\n"
done >> demulti_dna/demultiplex_summary.tsv
ls $PWD/demulti_dna/demultiplex_summary.tsv
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
for DataDir in $(ls -d demulti_dna/*/*/*/*/* | tail -n+2); do
  for R1 in $(ls $DataDir/*_R1*.fq | grep -v 'ambiguous'); do
    R2=$(echo $R1 | sed 's/_R1/_R2/g')
    echo $(basename $R1)
    echo $(basename $R2)
    Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    while [ $Jobs -gt 10 ]; do
      sleep 10s
      printf "."
      Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
    done
    printf "\n"
    Locus=$(echo $R1 | rev | cut -f1 -d '_' | rev | tr -d '.fq')
    Prefix=$(echo $DataDir | cut -f3,4,5,6 -d '/' | sed 's&/&_&g')"_${Locus}"
    OutDir="processed_dna/"$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
    # echo $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/2020_analysis
    qsub $ProgDir/sub_process_reads_2020.sh $R1 $R2 $Locus $Prefix $OutDir
  done
done
```


Summarise reads merged:
```bash
printf "RunName\tPlate\tLocus\tSoil\tDilution\tRep\t16S\tITS\tTef\tSIX5\tT4\tT6\n" > processed_dna/merged_summary.tsv
for RunDir in $(ls -d demulti_dna/*/*/*/*/*); do
  Plate=$(echo $RunDir | rev | cut -f5 -d '/' | rev)
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Soil=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName="processed_dna/$Plate/$Locus/$Soil/$Dilution/merged/${Locus}_${Soil}_${Dilution}_${Rep}*"
  T0Lines=$(cat ${RunName}*_*16S_prefilter.fa | grep '>'| wc -l)
  ItsLines=$(cat ${RunName}*_*ITS_prefilter.fa | grep '>'| wc -l)
  TefLines=$(cat ${RunName}*_*TEF_prefilter.fa | grep '>'| wc -l)
  Six5Lines=$(cat ${RunName}*_*SIX5_prefilter.fa | grep '>'| wc -l)
  T4Lines=$(cat ${RunName}*_*OG13890_prefilter.fa | grep '>'| wc -l)
  T6Lines=$(cat ${RunName}*_*OG4952_prefilter.fa | grep '>'| wc -l)
  printf "${Locus}_${Soil}_${Dilution}_${Rep}\t$Plate\t$Locus\t$Soil\t$Dilution\t$Rep\t$T0Lines\t$ItsLines\t$TefLines\t$Six5Lines\t$T4Lines\t$T6Lines\n"
done >> processed_dna/merged_summary.tsv
ls $PWD/processed_dna/merged_summary.tsv
```


Summise reads filtered:
```bash
printf "RunName\tPlate\tLocus\tSoil\tDilution\tRep\t16S\tITS\tTef\tSIX5\tT4\tT6\n" > processed_dna/filtered_summary.tsv
for RunDir in $(ls -d demulti_dna/*/*/*/*/*); do
  Plate=$(echo $RunDir | rev | cut -f5 -d '/' | rev)
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Soil=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Dilution=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  Rep=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName="processed_dna/$Plate/$Locus/$Soil/$Dilution/filtered/${Locus}_${Soil}_${Dilution}_${Rep}*"
  # RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  T0Lines=$(cat ${RunName}*_*16S.filtered.fa | grep '>'| wc -l)
  ItsLines=$(cat ${RunName}*_*ITS.filtered.fa | grep '>'| wc -l)
  TefLines=$(cat ${RunName}*_*TEF.filtered.fa | grep '>'| wc -l)
  Six5Lines=$(cat ${RunName}*_*SIX5.filtered.fa | grep '>'| wc -l)
  T4Lines=$(cat ${RunName}*_*OG13890.filtered.fa | grep '>'| wc -l)
  T6Lines=$(cat ${RunName}*_*OG4952.filtered.fa | grep '>'| wc -l)
  printf "${Locus}_${Soil}_${Dilution}_${Rep}\t$Plate\t$Locus\t$Soil\t$Dilution\t$Rep\t$T0Lines\t$ItsLines\t$TefLines\t$Six5Lines\t$T4Lines\t$T6Lines\n"
done >> processed_dna/filtered_summary.tsv
ls $PWD/processed_dna/filtered_summary.tsv
```

## OTU assignment
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. vsearch is used for dereplication and also provides clustering outputs.

 * All read files associated with the each locus in the project are concatenated
 * Clustering run on all the data for this locus within the project
 * Quantification can then be performed for each treatment against the total set

```bash
  # Concatenate files
  for Locus in 16S ITS TEF SIX5 OG13890 OG4952; do
  # for Locus in 16S; do
  # for Locus in SIX5; do
    OutDir=clustering/$Locus
    mkdir -p $OutDir
    # cat processed_dna/*/*/*/*/filtered/*_${Locus}.filtered.fa >> $OutDir/${Locus}_concatenated.temp.fa
    ls -lh $OutDir/${Locus}_concatenated.temp.fa
    # Submit clustering
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts/2020_analysis
    qsub $ProgDir/submit_clustering_2020.sh $OutDir/${Locus}_concatenated.temp.fa $OutDir $Locus
  done
```


### Assign taxonomy

https://www.drive5.com/usearch/manual/utax_or_sintax.html
Taxonomy can be assigned using SINTAX or UTAX methods:

```
UTAX and SINTAX have different strengths and weaknesses. On short 16S tags like V4, SINTAX and RDP have very similar performance. On longer 16S sequences and on ITS sequences, SINTAX is better than RDP. SINTAX is simpler because it doesn't need training, while training UTAX or RDP is quite challenging if you want to use your own database. UTAX is the only algorithm which tries to account for sparse reference data and has the lowest over-classification rate of any algorithm (except possibly the k-nearest-neighbor method in mothur, but knn has low sensitivity in general). However, UTAX sometimes has lower sensitivity than SINTAX to known taxa. Neither algorithm is a clear winner over the other.
```

* Greg uses Sintax as it is A) recomended by the author B) easier to make custom databases with
* An alternative to database-based taxonomy assignment is to used BLAST. It may be worth investigating automated BLAST lookups vs NCBI.

usearch hosts databases for some loci. Silva also has a database for the stremenophiles (oomycetes) that can be downloaded. Otherwise, personal databases must be made.
https://unite.ut.ee/repository.php
https://www.arb-silva.de/browser/

Reference databases were made for each locus in the project directory:

The 16S, ITS, TEF, OG13890 and OG4952 databases were already installed on the cluster.
These were copied into the project directory.

These databases were set up as part of previous ampseq runs (see 02_artificaial_pools_run2.md and 03_soil_samples_plate3.md)

```bash
cp -r /home/groups/harrisonlab/project_files/fusarium_ampseq/databases .
```

A SIX5 database was also created from a fasta file curated on my local computer in GENEIOUS.

Important - It was found that although the PCR primers contained the sequence CTGATGGCAAAGGTCATAGAATGTT, SIX5 reads once sequenced and processed were missing
the expected terminal 'A' as if the primer sequence CTGATGGCAAAGGTCATAGAATGTTT had
been trimmed off of the sequence (it wasn't). Could this be a PCR artifact? The
SIX5 fasta database was edited to contain the expected and edited sequence accordingly.

Important - When the two SIX5 amplicons were used (or the two expected amplicons
and two observed amplicons - discussed above), taxonomy assignment failed for OTUs
vs the database despite sequences sharing 100% identity. This was investigated
and wasn't found to be a result of newline characters etc. but a requirement to
have a greater number of sequences that show some homology to the expected
organism. For example, combining the SIX5 target sequences with OG4952 sequences
into a single database did not aid the identification of the zOTU sequence as
the SIX5 target (despite 100% identity). However making a hybrid sequence of an
OG4952 sequence and adding it into the SIX5 targets (labelling it as a control
sequence) allowed correct identification of the zOTU as the SIX5 sequence that
it shared 100% homology with. I conclude that this is to do with bootstrapping
support for identification of sequences as taxa, and taxa that have very
distinct sequences from other entries in the database may cause problems. This
may also explain why Trichoderma was missing in the previous project's TEF analysis.
See: https://www.drive5.com/usearch/manual/sintax_algo.html

```bash
DbDir=databases/SIX5
mkdir -p $DbDir
# assembly from my lcoal computer was copied here using the command:
# scp /Volumes/GGB/Harrison/Projects/AHDB-FOC-2019-C300106/Science/Obj4-Establish-Potential-AmpliconSeq/SIX5_primer_design/SIX5_amplicon_nr_db.fasta cluster:/projects/fusarium_ampseq/databases/SIX5/.

cd $DbDir
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -makeudb_sintax SIX5_amplicon_nr_db.fasta -output SIX5_amplicon_nr_db.udp
cd -
```


```bash
# for Locus in 16S ITS TEF SIX5 OG13890 OG4952; do
# for Locus in 16S ITS TEF OG13890 OG4952; do
for Locus in 16S; do
# for Locus in TEF SIX5 OG13890 OG4952; do
  for Prog in usearch vsearch; do
    for Type in zOTUs; do
      OtuFa=$(ls clustering/$Locus/${Locus}_cluster_${Prog}_${Type}.fa)
      RefDb=$(ls databases/$Locus/*.udp)
      Prefix=$(basename ${OtuFa%.fa})
      OutDir=$(dirname $OtuFa)
      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
      # echo "$OtuFa $RefDb $Prefix $OutDir"
      qsub $ProgDir/submit_taxonomy.sh $OtuFa $RefDb $Prefix $OutDir
    done
  done
done
```


## Quantification

Quantification of infested field soil samples from the 2019 project

```bash
  # ITS mix (ITS and 16S)
  # for Locus in 16S; do
  for Locus in ITS 16S; do
    Threshold=114
    Identity=0.97
    Prefix="${Locus}_2019_samples"
    echo $Prefix
    OutDir=quantified/$Locus
    mkdir -p $OutDir
    for OtuType in zOTUs; do
      for RefDb in $(ls clustering/$Locus/*_cluster_usearch_zOTUs_taxa.fa); do
        QueryReads=$OutDir/${Prefix}_reads_appended.fa
        cat processed_dna/1/ITS_mix/FO*/NA/filtered/ITS_mix_FO*_NA_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_1_${Locus}./g" | cut -f1 -d '.' > $QueryReads
        cat processed_dna/2/ITS_mix/FO*/NA/filtered/ITS_mix_FO*_NA_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_2_${Locus}./g" | cut -f1 -d '.' >> $QueryReads
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
      done
    done
  done
  # TEF mix (TEF, SIX5, OG13890 and OG4952)
  for Locus in TEF SIX5 OG13890 OG4952; do
  # for Locus in SIX5 OG13890 OG4952; do
    Threshold=114
    Identity=0.97
    Prefix="${Locus}_2019_samples"
    echo $Prefix
    OutDir=quantified/$Locus
    mkdir -p $OutDir
    for OtuType in zOTUs; do
      for RefDb in $(ls clustering/$Locus/*_cluster_usearch_zOTUs_taxa.fa); do
        QueryReads=$OutDir/${Prefix}_reads_appended.fa
        cat processed_dna/1/TEF_mix/FO*/NA/filtered/TEF_mix_FO*_NA_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_1_${Locus}./g" | cut -f1 -d '.' > $QueryReads
        cat processed_dna/2/TEF_mix/FO*/NA/filtered/TEF_mix_FO*_NA_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_2_${Locus}./g" | cut -f1 -d '.' >> $QueryReads
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
      done
    done
  done
  # The diagnostic loci were also run requiring 100% read identity to a taxon
  for Locus in TEF SIX5 OG13890 OG4952; do
  # for Locus in SIX5 OG13890 OG4952; do
    Threshold=114
    Identity=1.0
    Prefix="${Locus}_2019_samples_no-mismatches"
    echo $Prefix
    OutDir=quantified/$Locus
    mkdir -p $OutDir
    for OtuType in zOTUs; do
      for RefDb in $(ls clustering/$Locus/*_cluster_usearch_zOTUs_taxa.fa); do
        QueryReads=$OutDir/${Prefix}_reads_appended.fa
        cat processed_dna/1/TEF_mix/FO*/NA/filtered/TEF_mix_FO*_NA_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_1_${Locus}./g" | cut -f1 -d '.' > $QueryReads
        cat processed_dna/2/TEF_mix/FO*/NA/filtered/TEF_mix_FO*_NA_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_2_${Locus}./g" | cut -f1 -d '.' >> $QueryReads
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
      done
    done
  done
```

Quantification of spiked soil samples as part of the 2020 project

```bash
# for Locus in 16S; do
for Locus in ITS 16S; do
  Threshold=114
  Identity=0.97
  Prefix="${Locus}_soil_samples"
  echo $Prefix
  OutDir=quantified/$Locus
  mkdir -p $OutDir
  for OtuType in zOTUs; do
    for RefDb in $(ls clustering/$Locus/*_cluster_usearch_zOTUs_taxa.fa); do
      QueryReads=$OutDir/${Prefix}_reads_appended.fa
      cat processed_dna/*/ITS_mix/S*/D*/filtered/ITS_mix_S*_D*_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_${Locus}./g" | cut -f1 -d '.' > $QueryReads
      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
      qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
    done
  done
done

  # TEF mix (TEF, SIX5, OG13890 and OG4952)
  # for Locus in TEF SIX5 OG13890 OG4952; do
  for Locus in TEF SIX5 OG13890 OG4952; do
  # for Locus in SIX5 OG13890 OG4952; do
    Threshold=114
    Identity=0.97
    Prefix="${Locus}_soil_samples"
    echo $Prefix
    OutDir=quantified/$Locus
    mkdir -p $OutDir
    for OtuType in zOTUs; do
      for RefDb in $(ls clustering/$Locus/*_cluster_usearch_zOTUs_taxa.fa); do
        QueryReads=$OutDir/${Prefix}_reads_appended.fa
        cat processed_dna/*/TEF_mix/S*/D*/filtered/TEF_mix_S*_D*_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_${Locus}./g" | cut -f1 -d '.' > $QueryReads
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
      done
    done
  done
  # The diagnostic loci were also run requiring 100% read identity to a taxon
  for Locus in TEF SIX5 OG13890 OG4952; do
    Threshold=114
    Identity=1.0
    Prefix="${Locus}_soil_samples_no-mismatches"
    echo $Prefix
    OutDir=quantified/$Locus
    mkdir -p $OutDir
    for OtuType in zOTUs; do
      for RefDb in $(ls clustering/$Locus/*_cluster_usearch_zOTUs_taxa.fa); do
        QueryReads=$OutDir/${Prefix}_reads_appended.fa
        cat processed_dna/*/TEF_mix/S*/D*/filtered/TEF_mix_S*_D*_*_${Locus}.filtered.fa  | sed "s/_${Locus}/_${Locus}./g" | cut -f1 -d '.' > $QueryReads
        ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
        qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir $Threshold $Identity
      done
    done
  done
```

## Generation of plots

Plots were generated using R-studio from my local computer using data downloaded
onto EMQA

### FOC FON FOM soils

Documented in:
* 2020_analysis/plot_FOC_FON_FOM_soils_ITS.r
