# Bioinformatic pipeline

## Installation
Install required softwares
- [QIIME2 amplicon version 2026.1](https://library.qiime2.org/quickstart/amplicon)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://seqera.io/multiqc/)
- [R](https://cran.rstudio.com/)
- [R studio](https://posit.co/download/rstudio-desktop/) (not mandatory but suggested)


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

> [!IMPORTANT]
> The default values of the dada2 parameters are usually good for most applications, but you should experiment with the values of `--p-trunc-len-f` and `--p-trunc-len-r` and compare the results (in terms of model fitting, number of retained sequences per sample and sequences length) to choose the best values. Here an example of how to do it.

```bash
# truncation values for forward and reverse reads to be tested
TruncLenF=(250 250 240 240 230 230 220)
TruncLenR=(250 240 240 230 230 220 220)

for i in "${!TruncLenF[@]}"
do

  TLF=${TruncLenF[i]}
  TLR=${TruncLenR[i]}
    
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs approccio_GREZZO/${DIR}/seqs.qza \
    --p-trim-left-f $FW_LEN \
    --p-trim-left-r $RV_LEN \
    --p-trunc-len-f $TLF \
    --p-trunc-len-r $TLR \
    --p-max-ee-f 2 \
    --p-max-ee-r 2 \
    --p-trunc-q 2 \
    --p-n-reads-learn 1000000 \
    --p-pooling-method 'pseudo' \
    --p-n-threads $JOBS \
    --o-table table_${TLF}_${TLR}.qza \
    --o-representative-sequences rep-seqs_${TLF}_${TLR}.qza \
    --o-denoising-stats denoising-stats_${TLF}_${TLR}.qza \
    --o-base-transition-stats base-transition-stats_${TLF}_${TLR}.qza \
    --verbose

  # visualize denoising stats
  qiime metadata tabulate \
    --m-input-file denoising-stats_${TLF}_${TLR}.qza \
    --o-visualization denoising-stats_${TLF}_${TLR}.qzv

  # visualize base transition stats
  qiime dada2 plot-base-transitions \
    --i-base-transition-stats base-transition-stats_${TLF}_${TLR}.qza \
    --o-visualization base-transition-stats_${TLF}_${TLR}.qzv

  # Export files for checking them in R
  tmpdir="tmp_${TLF}_${TLR}"
  mkdir -p "$tmpdir"
  
  qiime tools export --input-path "denoising-stats_${TLF}_${TLR}.qza" --output-path "$tmpdir/exported_stats"
  cp "${tmpdir}/exported_stats/stats.tsv" "denoising-stats_${TLF}_${TLR}.tsv"

  qiime tools export --input-path "table_${TLF}_${TLR}.qza" --output-path "$tmpdir/exported_table"
  biom convert -i "$tmpdir/exported_table/feature-table.biom" -o "table_${TLF}_${TLR}.tsv" --to-tsv

  qiime tools export --input-path "rep-seqs_${TLF}_${TLR}.qza" --output-path "$tmpdir/exported_seqs"
  cp "${tmpdir}/exported_seqs/dna-sequences.fasta" "rep-seqs_${TLF}_${TLR}.fasta"

  rm -r "$tmpdir"

done

```

In R you can use this code for checking the dada2 output.

```bash
### needed libraries (you'll need to install them the first time using 'install.packages(<"package_name">)')
library(tidyr)
library(dplyr)
library(ggpubr)
library(Biostrings)
library(plotly)

### PLOT 1 - for checking number of retained reads

# read data
cc=0
dat_list <- list()
for (fil in list.files(pattern = "denoising-stats_.*tsv")) {
  cc=cc+1
  xx<-read.csv(fil, sep="\t")
  xx<-xx[-1,]
  xx_long <- xx %>% 
    pivot_longer(
      cols      = -sample.id,
      names_to  = "metric",
      values_to = "value"
    )
  dat_list[[cc]]<-as.data.frame(xx_long)
}
names(dat_list) <- sub("^[^_]+_([^\\.]+)\\.tsv$", "\\1",
                       list.files(pattern = "denoising-stats_.*tsv"))

df <- bind_rows(dat_list, .id = "source_df")
df$value<-as.numeric(df$value)

# define comparisons anc columns of interest
my_comparisons<-combn(unique(df$source_df), 2, simplify = FALSE)

ToPlot<-c("percentage.of.input.passed.filter",
          "percentage.of.input.merged",
          "percentage.of.input.non.chimeric")

# make the plot
plts<-list()
for (i in 1:length(ToPlot)) {
  dd<-df[df$metric==ToPlot[i],]
  plts[[i]]<-ggboxplot(dd, x = "source_df", y = "value",
                        add = "jitter",
                        title = ToPlot[i]) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggarrange(plotlist = plts, ncol=3)
#ggsave("Figure_xx.png", device="png", 
#       width = 29.7*1.5, height = 21, units = "cm", bg="white")

### PLOT 2 - for checking amplicon length

# read data
SeqsLengthTotalData<-list()
SeqsLengthPerSample<-list()
for (fil in list.files(pattern="^table_[0-9_]+\\.tsv$")) {
  id <- sub("^table_(.*)\\.tsv$", "\\1", fil)
  
  # read table
  counts <- read.delim(fil, skip = 1, check.names = FALSE) 
  rownames(counts) <- counts$"#OTU ID"
  counts <- counts[ , -1, drop = FALSE]
  
  # read sequences
  seqs <- readDNAStringSet(paste0("rep-seqs_",id,".fasta"))
  seq_lengths <- width(seqs)
  names(seq_lengths) <- names(seqs)
  
  # match only shared features (just for safety)
  shared_features <- intersect(rownames(counts), names(seq_lengths))
  counts <- counts[shared_features, , drop = FALSE]
  seq_lengths <- seq_lengths[shared_features]
  
  # save length of each seqs in total data
  seq_lengths2<-as.data.frame(seq_lengths)
  seq_lengths2$ASV<-row.names(seq_lengths2)
  SeqsLengthTotalData[[id]]<-seq_lengths2
  
  # compute average length per sample
  average_lengths <- sapply(counts, function(sample_counts) {
    total_reads <- sum(sample_counts)
    if (total_reads == 0) return(NA)
    sum(sample_counts * seq_lengths) / total_reads
  })
  
  # save average length per sample
  SeqsLengthPerSample[[id]]<-as.data.frame(average_lengths)
  SeqsLengthPerSample[[id]]$SampleID<-row.names(SeqsLengthPerSample[[id]])
  
}

# make the plot - static version
PLOTS_length<-list()
for (i in 1:length(SeqsLengthTotalData)) {
  df<-SeqsLengthTotalData[[i]]
  PLOTS_length[[i]] <- ggplot(df, aes(x = seq_lengths)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    labs(title = names(SeqsLengthTotalData)[i],
         x = "Sequence Length",
         y = "Density") +
    theme_pubr()
}

ggarrange(plotlist = PLOTS_length, nrow = 3, ncol = 3)
#ggsave("Figure_xx.png", device="png", 
#       width = 29.7*1.5, height = 21, units = "cm", bg="white")

# make the plot - interactive version
long_df <- bind_rows(
  lapply(names(SeqsLengthTotalData), function(name) {
    df <- SeqsLengthTotalData[[name]]
    df$method <- name
    return(df)
  })
)

p<-ggplot(long_df, aes(x = seq_lengths, color = method)) +
  geom_density(size = 1) +
  theme_minimal() +
  labs(
    x = "Sequence Length",
    y = "Density",
    color = "Denoising Method",
    title = "Sequence Length Distribution by Denoising Method"
  )

ggplotly(p)

```

Here some examples of the files produced:
- [denoising stats](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/denoising-stats.qzv)
- [base transition stats](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/base-transition-stats.qzv)
- [PLOT1](https://leafbeetle.github.io/meta_pipeline/examples/PLOT1.png)
- [PLOT2 - static version](https://leafbeetle.github.io/meta_pipeline/examples/PLOT2.png)
- [PLOT2 - interactive version](https://leafbeetle.github.io/meta_pipeline/examples/PLOT2_interactive.html)


After selecting the best dada2 results chenge their name for simplicity.
```bash
TLF=<selcted_TLF>
TLR=<selected_TLR>

cp rep-seqs_${TLF}_${TLR}.qza rep_seqs.qza
cp table_${TLF}_${TLR}.qza table.qza 

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
Here examples of these visualizations:
- [table](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/table.qzv)
- [sequences](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/rep-seqs.qzv)


## Taxonomic classification
You can download the silva database for the bacterial 16S [here](https://drive.google.com/drive/folders/1eNXLXgK7321b_EBJExokK6vjuLgrRAf2?usp=sharing).

To optimize taxonomic classification and reduce database complexity we can trim the reference sequences to contain only the region actually amplified by the primers used in this study.
```bash
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer <sequence_of_forward_primer> \
  --p-r-primer <sequence_of_reverse_primer> \
  --p-identity 0.8 \
  --p-min-length 200 \
  --p-max-length 550 \
  --p-n-jobs $JOBS \
  --p-read-orientation forward \
  --o-reads silva-138-99-seqs_<primer_pair_name>.qza

qiime rescript dereplicate \
  --i-sequences silva-138-99-seqs_<primer_pair_name>.qza \
  --i-taxa silva-138-99-tax.qza \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138-99-seqs_<primer_pair_name>_derep.qza \
  --o-dereplicated-taxa silva-138-99-tax_<primer_pair_name>_derep.qza \
  --p-threads $JOBS
```

Now we are ready for the taxonomic classification. First of all train a classifier on the SILVA reference database that we built before.
```bash
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs_<primer_pair_name>_derep.qza \
  --i-reference-taxonomy silva-138-99-tax_<primer_pair_name>_derep.qza \
  --p-verbose \
  --o-classifier NBclassifier.qza
```

Now we can use the trained classifier to assign a taxonomic label to each ASV. The `--p-confidence` parameter is used for setting a confidence threshold for limiting taxonomic depth, in this case we set it to 70% that is usually a good value for bacterial 16S.
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier NBclassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.70 \
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
[Here](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/barplot.qzv) an example of a barplot visualization.













