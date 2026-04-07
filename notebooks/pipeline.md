# Bioinformatic pipeline

## Installation
Install required softwares
- [QIIME2 amplicon version 2026.1](https://library.qiime2.org/quickstart/amplicon)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://seqera.io/multiqc/)


## Preliminary settings

We start defining some variables so that the pipeline can be applied to different datasets.
- `RAWDIR` for the directory containing raw sequences (i.e. the fastq files)
- `WORKDIR` for the working directory where we want to save the results of the analyses
- `JOBS` for the maximum number of processes that can be run simultaneously (i.e. maximum number of available cores)
- `ENV` for the name of the conda environment where we installed QIIME2
- `FW_LEN` for the length of the forward primer
- `RV_LEN` for the length of the reverse primer

```bash
RAWDIR=<path/to/rawdata/dir>

WORKDIR=<path/to/outputs>

JOBS=<number_of_cores>

ENV=qiime2-amplicon-2026.1

FW_LEN=<length_of_forward_primer>
RV_LEN=<length_of_reverse_primer>

```


## Reads quality check
First of all check the quality of the raw data using the softwares **FastQC** and **MultiQC**. 

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
Look at the multiqc report for checking sequence quality. Here an [example](https://MontagnaLab.github.io/InnovativeApproachesForInvertebrateBiodiversityMonitoring/data/multiqc_report.html).


## Import sequences in QIIME2
Make a manifest file named `manifest.tsv` using the SampleID of the metadata file and save it in the working directory `WORKDIR`. 
It should be a TSV (tab separated values) file with three columns: sample ID, path to forward reads, path to reverse reads.
Like these examples.

- paired-end reads (e.g., Illumina)
```bash
sample-id    forward-absolute-filepath    reverse-absolute-filepath
S0055    <path/to>/SRR10896373_AF33_1.fastq.gz	<path/to>/SRR10896373_AF33_2.fastq.gz
S0056    <path/to>/SRR10896374_AF32_1.fastq.gz	<path/to>/SRR10896374_AF32_2.fastq.gz
...
```

- single-end reads (e.g., IonTorrent)
```bash
sample-id     absolute-filepath
S0001	<path/to>/SRR13222562_gut_microbiota_of_Leptinotarsa_decemlineata_from_China_Xinjiang_Wenquan_County.fastq.gz
S0002	<path/to>/SRR13222563_gut_microbiota_of_Leptinotarsa_decemlineata_from_China_Xinjiang_Wenquan_County.fastq.gz
```

Import sequences using the manifest file.
```bash
cd $WORKDIR

# for paired-end reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path manifest.tsv \
  --output-path seqs.qza

# for single-end reads
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path manifest.tsv \
  --output-path seqs.qza
```

To see a summary of the imported sequences we can create a visualization.
```bash
qiime demux summarize \
  --i-data seqs.qza \
  --o-visualization seqs.qzv
```
> [!TIP]
> QIIME2 visualizations can be visualized using the command `qiime tools view <filename.qzv>` or loading the files in the [online visualizer](https://view.qiime2.org/).

Look at this visualization for an overall report on sequence quality. Here an [example](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/seqs.qzv).


## Denoising

We are now ready for removing non-biological variation from our data. We use the [DADA2 algorithm](https://benjjneb.github.io/dada2/) implemented in [q2-dada2](https://github.com/qiime2/q2-dada2) that models and corrects sequencing errors to infer exact biological sequences (amplicon sequence variants, ASVs). The most important parameters are:
- `--p-trim-left-f` and `--p-trim-left-r`, corresponding to the length of the forward and reverse primer respectively
- `--p-trunc-len-f` and `--p-trunc-len-r`, corresponding to the length at which to trunc the sequences due to quality drop
- `--p-max-ee-f` and `--p-max-ee-r`, reads with number of expected errors higher than this value will be discarded
- `--p-trunc-q`, reads are truncated at the first instance of a quality score less than or equal to this value
- `--p-n-reads-learn`, number of reads used for training the error model, 1M is usually enough, but it may be increased for big datasets
```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs seqs_trimmed.qza \
  --p-trim-left-f <length_of_forward_primer> \
  --p-trim-left-r 18 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 200 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-n-reads-learn 1000000 \
  --p-pooling-method 'pseudo' \
  --p-n-threads $JOBS \
  --o-table table.qza \
  --o-representative-sequences rep-seqs_MixedOrientation.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza \
  --verbose

# visualize denoising stats
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# visualize base transition stats
qiime dada2 plot-base-transitions \
  --i-base-transition-stats base-transition-stats.qza \
  --o-visualization base-transition-stats.qzv
```

> [!IMPORTANT]
> You should experiment with the values of `--p-trunc-len-f` and `--p-trunc-len-r` parameters and compare the results (in terms of number of retained sequences per sample and sequences length) to choose the best values.

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














################################################################################################################################################à

#### 2.3.4. Extract the amplified region from the reference database

To optimize taxonomic classification and reduce database complexity we can **trim the reference sequences**
to contain only the region actually amplified by the primers used in this study.
```bash
qiime feature-classifier extract-reads \
  --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
  --p-f-primer ASCYGYGGTAAYWCCAGC \
  --p-r-primer TCHNHGNATTTCACCNCT \
  --p-identity 0.8 \
  --p-min-length 200 \
  --p-max-length 550 \
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













