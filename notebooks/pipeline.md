# Bioinformatic pipeline

## Install QIIME2
Install version 2026.1 of QIIME2 amplicon following the instructions [here](https://library.qiime2.org/quickstart/amplicon).

## Preliminary settings

We start defining some variables so that the pipeline can be applied to different datasets.
- `RAWDIR` for the directory containing raw sequences (i.e. the fastq files)
- `WORKDIR` for the working directory where we want to save the results of the analyses
- `JOBS` for the maximum number of processes that can be run simultaneously (i.e. maximum number of available cores)
- 
```bash
RAWDIR=<path/to/rawdata/dir>

WORKDIR=<path/to/outputs>

JOBS=<number_of_cores>
```

Finally we define a variable containing the name of the conda environment where we installed **QIIME2**.
```bash
ENV=qiime2-amplicon-2026.1
```


### 2.2. Reads quality check
First of all we check the quality of the raw data using the softwares **FastQC** and **MultiQC**. 
If you installed these softwares you can run the commands, but since this can be a quite long step 
it may be better to look directly at the final report [here](https://MontagnaLab.github.io/InnovativeApproachesForInvertebrateBiodiversityMonitoring/data/multiqc_report.html).

```bash
# move to the directory where we want to save the outputs #
cd "$WORKDIR"
# create a directory for FastQC outputs #
mkdir -p FastQC_output

# run fastQC #
for file in "$RAWDIR"/*.fq.gz; do
    SAMPLE=$(basename "$file")
    fastqc -t $JOBS -o FastQC_output "$file"
    echo "Processed $SAMPLE"
done

# run multiqc #
cd FastQC_output
multiqc .
```

**FastQC** produces one report for each fastq file, **MultiQC** combines those reports in a single one for all samples.
Let's have a look at the reports.

⚠️
There are still a few illumina universal adapters in some reads, we will remove them later.

### 2.3 Create reference database
To taxonomically classify our sequences we need to construct a reference database.
In this case we will use [SILVA v.138.1](https://www.arb-silva.de/), a quality checked and regularly updated 
database of aligned small (16S/18S, SSU) and large subunit (23S/28S, LSU) ribosomal RNA (rRNA) sequences 
for all three domains of life (Bacteria, Archaea and Eukarya). 
This can be a very long step since the database is quite huge, so it may be better if you download the 
already prepared [reference sequences](https://github.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/blob/main/data/silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza) and [taxonomy](https://github.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/blob/main/data/silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza).

#### 2.3.1 Obtain and clean the reference sequences

First of all we activate the conda environment containing **QIIME2**.

```bash
conda activate $ENV
```

We will use the [q2-RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) plugin that has several built-in functions
for managing and curating reference sequence databases.

**Download** the SILVA RNA sequences and the associated taxonomic labels.
```bash
qiime rescript get-silva-data \
  --p-version '138.1' \
  --p-target 'SSURef_NR99' \
  --p-include-species-labels \
  --p-no-rank-propagation \
  --parallel \
  --o-silva-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
  --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza \
  --verbose
```

**Reverse-transcribe** the RNA sequences to obtain the corresponding DNA sequences.
```bash
qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza
```

**Remove low quality sequences**, in this case those with 5 or more degenerated bases 
and/or containing homopolymers with 8 or more bases)
```bash
qiime rescript cull-seqs \
  --i-sequences silva-138.1-ssu-nr99-seqs.qza \
  --p-num-degenerates 5 \
  --p-homopolymer-length 8 \
  --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
  --p-n-jobs $JOBS # number of concurrent processes
```

#### 2.3.2. Clean the reference taxonomy
Even if SILVA is a curated database the taxonomic labels can be a bit messy, we will use **R** to clean them a bit.

First we need to **export the taxonomy** in the TSV format (Tab Separated Values), so that **R** can read it.
```bash
# export taxonomy to tsv file
qiime tools export \
  --input-path silva-138.1-ssu-nr99-tax-cleaned.qza \
  --output-path exp
```

Now **open Rstudio or an R terminal** and set the working directory to the same of `$WORKDIR`.

```r
setwd("<path/to/outputs>")
```

**Load the required libraries**, install them if they are not already installed. 
```r
#install.packages("dplyr")
library(dplyr)
#install.packages("stringr")
library(stringr)
```

**Read the file** containing the SILVA taxonomy.
```r
data <- read.table("exp/taxonomy.tsv", header=T, sep="\t")
```

**Remove bad species-level ID** that do not correspond to taxon names (i.e. those starting with lowercase) 
and fill empty genus rank using species-level identifications.
```r
CleanData <- data %>%
  mutate(Taxon = case_when(
      # Case 1: First letter after "s__" is lowercase
      str_detect(Taxon, "s__[a-z]") ~ str_replace(Taxon, "s__[a-z].*", "s__"),
      # Case 2: First letter after "s__" is uppercase
      str_detect(Taxon, "s__[A-Z]") ~ str_replace(Taxon, 
                                                  "g__; s__([^_]+)",  # Capture "g__;" and genus portion after "s__"
                                                  "g__\\1; s__\\1"),  # Append genus after "g__" and keep it after "s__"
      # Default: Leave other rows unchanged
      TRUE ~ Taxon
    )
  )
```

**Remove bad genus-level ID** that do not correspond to taxon names (i.e. those starting with lowercase).
```r
CleanData2 <- CleanData %>%
  mutate(
    Taxon = case_when(
      # If the first letter after "g__" is lowercase, remove everything after "g__" and keep "s__"
      str_detect(Taxon, "g__[a-z]") ~ str_replace(Taxon, "g__[a-z].*s__", "g__; s__"),
      
      # Otherwise, leave the string as it is
      TRUE ~ Taxon
    )
  )
```

**Save** the cleaned taxonomy to a TSV file.
```r
write.table(CleanData2, "silva-138.1-ssu-nr99-tax-cleaned.tsv", sep="\t", quote = F, row.names = F)
```

Close Rstudio or the R terminal and **go back to the bash terminal** with the QIIME2 environment
to import the cleaned SILVA taxonomy as a QIIME2 artifact.
```bash
qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-path silva-138.1-ssu-nr99-tax-cleaned.tsv \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path silva-138.1-ssu-nr99-tax-cleaned.qza
```


#### 2.3.3. Filter and dereplicate the reference database

Next we **remove** database entries with **too short SSU sequences** depending on the taxon.
In this case we use a minimum length threshold of 900 bp for Archaea, 1200 bp for Bacteria and 1400 bp for Eukaryota.
```bash
qiime rescript filter-seqs-length-by-taxon \
  --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
  --i-taxonomy silva-138.1-ssu-nr99-tax-cleaned.qza \
  --p-labels Archaea Bacteria Eukaryota \
  --p-min-lens 900 1200 1400 \
  --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
  --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza
```

Next, we **dereplicate entries** with the same taxonomic ID to remove redundance.
```bash
qiime rescript dereplicate \
  --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza \
  --i-taxa silva-138.1-ssu-nr99-tax-cleaned.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
  --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
  --p-threads $JOBS
```


#### 2.3.4. Extract the amplified region from the reference database

To optimize taxonomic classification and reduce database complexity we can **trim the reference sequences**
to contain only the region actually amplified by the primers used in this study.
```bash
qiime feature-classifier extract-reads \
  --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
  --p-f-primer ASCYGYGGTAAYWCCAGC \
  --p-r-primer TCHNHGNATTTCACCNCT \
  --p-identity 0.8 \
  --p-min-length 150 \
  --p-max-length 450 \
  --p-n-jobs $JOBS \
  --p-read-orientation forward \
  --o-reads silva-138.1-ssu-nr99-seqs_Euk575-895.qza
```

**Dereplicate** again (since some of these shorter sequences may be identical now) to remove redundance.
```bash
qiime rescript dereplicate \
  --i-sequences silva-138.1-ssu-nr99-seqs_Euk575-895.qza \
  --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --o-dereplicated-taxa silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --p-threads $JOBS
```


### 2.4. Import sequences in QIIME2
We are now ready to import the LUCAS sequences in QIIME2.

#### 2.4.1. Make a manifest file containing paths to raw reads
We need to create a manifest file telling QIIME2 where to find forward and reverse reads for each sample.
This is a TSV file with three columns: sample ID, path to forward reads, path to reverse reads.
First move to the directory containing the fastq files.
```bash
cd "$RAWDIR"
```

Get the **absolute paths to forward and reverse reads** and save it in separate files.
```bash
ls -1 "$PWD/"*1.fq.gz > R1.txt # paths to forward reads
ls -1 "$PWD/"*2.fq.gz > R2.txt # paths to reverse reads
```

Next we are going to extract **sample names** from the names of the fastq files.
The file are named following this scheme: "Lucas\<sampleID\>.\<sequencing-direction\>.fq.gz" 
(e.g., Lucas0001.1.fq.gz, Lucas0001.2.fq.gz, Lucas0002.1.fq.gz, Lucas0002.2.fq.gz, ...).
So the first 9 characters of the file names correspond to the sample names. 
Let's use some bash commands to extract sample names. 
```bash
# from each fastq file name extract the first 9 characters, sort them, and keep only unique values #
ls *fq.gz | cut -c 1-9 | sort | uniq > ids.txt
```

Combine sample names and absolute paths to make a manifest file.
```bash
# add column names to the manifest file #
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv

# add sample names and paths to the manifest file #
paste ids.txt R1.txt R2.txt >> manifest.tsv

# remove temporary files #
rm *.txt

# move the manifest file to the working directory #
mv manifest.tsv "$WORKDIR"/manifest.tsv
```

Let's have a look at the content of the manifest file to check that everything is ok.

#### 2.4.2. Import sequences
First move back to the working directory.
```bash
cd $WORKDIR
```

Now we can use the manifest file to import sequences in QIIME2 with `qiime tools import`. The main parameters of this command are:
- `--type` specifies the type of qiime2 artifact (QZA) to be created, you can use `qiime tools list-types` to see all the importable types available.
- `--input-format` specifies the format of the data to be imported, you can use `qiime tools list-formats` to see all the input data format available
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path manifest.tsv \
  --output-path seqs.qza 
```

To see a summary of the imported sequences we can create a visualization (QZV) with `qiime demux summarize`.
```bash
qiime demux summarize \
  --i-data seqs.qza \
  --o-visualization seqs.qzv
```
> [!TIP]
> QIIME2 QZV files can be visualized using the command `qiime tools view <filename.qzv>` or loading the files in the [online visualizer](https://view.qiime2.org/).

Let's have a look at [seqs.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/seqs.qzv).


If you remember the quality check reports there are still some Illumina adapter in the sequences, it is better to remove them using [q2-cutadapt](https://github.com/qiime2/q2-cutadapt). The command `qiime cutadapt trim-paired` has many parameters, in this case we just need `--p-adapter-f` and `--p-adapter-r` to specify the adapter sequence that we want to remove from the 3' end of forward and reverse reads, respectively.
```bash
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences seqs.qza \
  --p-cores $JOBS \
  --p-adapter-f AGATCGGAAGAG \
  --p-adapter-r AGATCGGAAGAG \
  --o-trimmed-sequences seqs_trimmed.qza \
  --verbose
```

### 2.5. Denoising and sequence re-orientation

We are now ready for removing non-biological variation from our data. We use the [DADA2 algorithm](https://benjjneb.github.io/dada2/) implemented in [q2-dada2](https://github.com/qiime2/q2-dada2) that models and corrects sequencing errors to infer exact biological sequences (amplicon sequence variants, ASVs). The most important parameters are:
- `--p-trim-left-f` and `--p-trim-left-r`, corresponding to the length of the forward and reverse primer respectively
- `--p-trunc-len-f` and `--p-trunc-len-r`, corresponding to the length at which to trunc the sequences due to quality drop
- `--p-max-ee-f` and `--p-max-ee-r`, reads with number of expected errors higher than this value will be discarded
- `--p-trunc-q`, reads are truncated at the first instance of a quality score less than or equal to this value
- `--p-n-reads-learn`, number of reads used for training the error model, 1M is usually enough, but it may be increased for big datasets
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs seqs_trimmed.qza \
  --p-trim-left-f 18 \
  --p-trim-left-r 18 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 200 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-n-reads-learn 1000000 \
  --p-pooling-method 'pseudo' \
  --p-n-threads $JOBS \
  --o-table table_MixedOrientation.qza \
  --o-representative-sequences rep-seqs_MixedOrientation.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza \
  --verbose
```

> [!IMPORTANT]
> In real life scenarios you should experiment with `--p-trunc-len-f` and `--p-trunc-len-r` parameters and compare the results (in terms of number of retained sequences per sample and sequences length) to choose the best values.

Now let's create visualizations for the two stats files produced by the DADA2 algorithm.
```bash
# visualize denoising stats
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# visualize base transition stats
qiime dada2 plot-base-transitions \
  --i-base-transition-stats base-transition-stats.qza \
  --o-visualization base-transition-stats.qzv
```
Let's have a look at [denoising-stats.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/denoising-stats.qzv) and [base-transition-stats.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/base-transition-stats.qzv).

Since in this study barcodes and adapters were added after PCR amplification each fastq file contained both forward and reverse reads. So sequences needs to be re-orientered using the reference database as guide. We can use the `qiime rescript orient-seqs` command from the [q2-RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) plugin for doing it.
```bash
qiime rescript orient-seqs \
  --i-sequences rep-seqs_MixedOrientation.qza \
  --i-reference-sequences silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --o-oriented-seqs rep-seqs.qza \
  --o-unmatched-seqs orientation_unmatched_sequences.qza
```

Now we can exclude unmatched sequences from the ASV table.
```bash
qiime feature-table filter-features \
  --i-table table_MixedOrientation.qza \
  --m-metadata-file orientation_unmatched_sequences.qza \
  --p-exclude-ids \
  --o-filtered-table table.qza
```

Create visualizations for the ASV table and ASV sequences files.
```bash
# table visualization
qiime feature-table summarize \
  --i-table table.qza \
  --m-metadata-file metadata.tsv \
  --o-feature-frequencies feature-frequencies.qza \
  --o-sample-frequencies sample-frequencies.qza \
  --o-summary table.qzv

# sequences visualization
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```
Let's have a look at [table.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/table.qzv) and [rep-seqs.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/rep-seqs.qzv).


### 2.6. Taxonomic classification
Now we are ready for the taxonomic classification. There are different methods in the [q2-feature-classfifier](https://github.com/qiime2/q2-feature-classifier?tab=readme-ov-file) plugin for this, this time we will use a machine learning approach with a naive bayes classifier. First of all let's train a classifier on the SILVA reference database that we built before.
```bash
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138.1-ssu-nr99-seqs_Euk575-895_derep-uniq.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax_Euk575-895_derep-uniq.qza \
  --p-verbose \
  --o-classifier NBclassifier.qza
```
> [!TIP]
> In real case scenarios, to improve taxonomic classification accuracy you may want to use also the `--i-class-weight` parameter to incorporate environment-specific taxonomic abundance information that can be generated by the [q2-clawback](https://github.com/BenKaehler/q2-clawback) plugin.

Now we can use the trained classifier to assign a taxonomic label to each ASV. The `--p-confidence` parameter is used for setting a confidence threshold for limiting taxonomic depth, in this case we set it to 90% that is usually a good value for metazoan.
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier NBclassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.90 \
  --p-n-jobs $JOBS \
  --o-classification taxonomy.qza
```

To visualize the taxonomic classification we can generate a barplot.
```bash
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization barplot.qzv
```
Let's have a look at [barplot.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/barplot.qzv).


### 2.7. Final filtering

As you can see in the barplot there are several non-invertebrate taxa, this is due to the use of primers that are not specific for invertebrates. Now let's filter the table to keep only the invertebrate taxa of our interest: nematodes, arthropods, tardigrads, annelids, rotifers, and flatworms. We can specify the taxa we want to keep using the `--p-include` parameter.
```bash
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include Nematozoa,Arthropoda,Tardigrada,Annelida,Rotifera,Platyhelminthes \
  --o-filtered-table invertebrates_table.qza
```

To reduce the noise we remove all the ASVs with less than 10 observations.
```bash
qiime feature-table filter-features \
  --i-table invertebrates_table.qza \
  --p-min-frequency 10 \
  --o-filtered-table invertebrates_table_clean.qza
```

Then generate a new barplot.
```bash
qiime taxa barplot \
  --i-table invertebrates_table_clean.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization invertebrates_barplot.qzv
```
Let's have a look at [invertebrates_barplot.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/invertebrates_barplot.qzv).













