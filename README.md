# For Job Application Review Purpose

Majority of the scripts are research working copies. They show the tools, parameters, and analysis choices used on real datasets.

## Contents

| Folder | Analysis Workflow | Related Publication |
| --- | --- | --- |
| [`16S-project`](16S-project) | QIIME 2 processing, followed by diversity and taxonomic analysis in R | **Title:** [Is Intestinal Microbiota Fully Restored After Chickens Have Recovered from Coccidiosis?](https://www.mdpi.com/2076-0817/14/1/81)<br>**Journal:** *Pathogens* (2025)<br>**Authors:** **Jiaqing Guo**, Zijun Zhao, Chace Broadwater, Isabel Tobin, Jing Liu, Melanie Whitmore, Guolong Zhang |
| [`nanopore-project`](nanopore-project) | Nanopore basecalling through taxonomic classification, followed by microbial community analysis in R | Manuscript in preparation |
| [`16S-database`](16S-database) | In-house two-step taxonomic classification using a custom reference database followed by RDP classification | — |
| [`metagenomic`](metagenomic) | Quality control, host read removal, assembly, gene prediction, and metagenomic binning | **Title:** [Dietary regulation on gut resistome linked with microbial amino acid metabolism in pigs](https://link.springer.com/article/10.1186/s42523-026-00542-z)<br>**Journal:** *Animal Microbiome* (2026)<br>**Authors:** Ziyu Liu, Chengming Meng, Junhua Shen, Huisong Wang, **Jiaqing Guo**, Jiangchao Zhao, Chunlong Mu, Weiyun Zhu |
| [`RNA-seq`](RNA-seq) | RNA-seq alignment, gene quantification, differential expression analysis, and GO and KEGG pathway enrichment analysis in R | **Title:** [Deoxycholic Acid Mitigates Necrotic Enteritis Through Selective Inhibition of Pathobionts and Enrichment of Specific Lactic Acid Bacteria](https://www.mdpi.com/2076-0817/14/7/688)<br>**Journal:** *Pathogens* (2025)<br>**Authors:** Melanie A. Whitmore, **Jiaqing Guo**, Dohyung M. Kim, Jing Liu, Isabel Tobin, Guolong Zhang |
| [`pythonNotes`](pythonNotes) | Python, pandas, and scikit-learn practice notebooks from coursework | — |

## 16S amplicon — `16S-project`

Paired-end 16S reads, processed in QIIME 2 and summarized in R. The biological comparison is treatment across a time course (days post infection).

**`16S_EMtimecourse_summer.sh`** — QIIME 2 (amplicon 2023.9 / 2025.4)

- Pull SRA runs and build a sample manifest
- Import reads, trim primers with cutadapt, merge pairs with vsearch
- Quality-filter, denoise 16S with Deblur, and build a MAFFT–FastTree phylogeny
- Assign taxonomy with a scikit-learn classifier, including a naive Bayes classifier trained on an RDP reference

**`16S_EMtimecourse_summer.R`** — phyloseq, vegan, tidyverse

- Clean the feature and taxonomy tables and drop rare features
- Alpha diversity (Observed, Shannon) with normality and variance checks
- Beta diversity: Bray–Curtis, Jaccard, UniFrac, and weighted UniFrac
- PCoA and PERMANOVA (`adonis2`) within each time point
- Relative-abundance bar plots by taxonomic rank

## Nanopore full-length 16S — `nanopore-project`

Full-length 16S from Oxford Nanopore, from raw signal to gut-segment community comparisons (crop through colon).

**`16S_fulength_misisibi_d42.sh`**

- Dorado basecalling (R10.4.1 sup model) and barcode demultiplexing
- Convert BAM to FASTQ; cutadapt to remove 16S primers and filter by length
- Taxonomy with EMU, including a project-specific reference built from NCBI accessions
- Optional clustering (CD-HIT), Medaka consensus polishing, and VSEARCH denoising

**`16s_fulength.R`** — vegan, tidyverse, ggplot2

- Fill incomplete taxonomy labels before aggregation
- Alpha diversity: Shannon, observed richness, Pielou evenness, with ANOVA or Kruskal–Wallis and compact letter displays
- Beta diversity: Bray–Curtis and Jaccard PCoA, overall and pairwise PERMANOVA
- Stacked taxonomy bar plots, rarefaction curves, and overlap (UpSet) across groups

## Custom 16S classifier — `16S-database`

**`run_classifier.sh`** is a small command-line wrapper around QIIME 2.

```bash
bash run_classifier.sh -w <work_dir> -d <database_dir> -c <rdp_classifier.qza>
```

- Classify representative sequences against a custom reference (vsearch consensus, 100% identity and query coverage)
- Keep hits labeled as EzBioCloud
- Send unassigned sequences to an RDP naive Bayes classifier
- Merge both sources into one taxonomy table and check that classification and search-hit tables have the same number of rows

## Shotgun metagenomics — `metagenomic`

**`Metagenomic.sh`** is the assembly and binning pipeline. Tools and versions are listed at the top of the script.

- Read QC with FastQC and MultiQC; trimming with fastp
- Remove host reads with Bowtie 2
- Assemble with MEGAHIT, both per sample and as a co-assembly
- Predict genes with Prodigal and cluster them with CD-HIT
- Search proteins with DIAMOND
- Bin contigs with MetaBAT2 from coverage profiles (CheckM2, MaxBin2, and CONCOCT are set up in the same environment)

## RNA-seq — `RNA-seq`

Chicken RNA-seq comparing Mock, DCA, NE, and DCA+NE.

**`RNA_Seq_Melanie_pipeline.sh`**

- FastQC / MultiQC, then fastp trimming
- HISAT2 index and alignment to the chicken genome
- Gene counts with featureCounts

**`RNA_Seq_Melanie.R`** — DESeq2, clusterProfiler, ggplot2

- Build a DESeq2 dataset, filter low-count genes, and fit `~ group`
- Contrasts: DCA vs Mock, NE vs Mock, DCA+NE vs Mock, DCA+NE vs NE
- Sample ordination, DEG counts, UpSet / Venn overlap
- Volcano plots and heatmaps of differentially expressed genes
- GO and KEGG enrichment, plus GSEA, with `org.Gg.eg.db` (chicken)

## Python notes — `pythonNotes`

**`Python_Notes.ipynb`** is a practice notebook (Coursera-style exercises), not a research analysis. It covers Python basics, pandas and NumPy, plotting, and a short scikit-learn section (including random forest). Small Python snippets inside the Nanopore shell script handle accession-to-taxonomy mapping and cluster bookkeeping.

## Tools

| Step | Tools |
| --- | --- |
| Amplicon | QIIME 2, cutadapt, vsearch, Deblur, scikit-learn / RDP classifiers |
| Long-read 16S | Dorado, samtools, cutadapt, EMU, CD-HIT, Medaka, VSEARCH |
| Metagenome | FastQC, fastp, MultiQC, Bowtie 2, MEGAHIT, Prodigal, CD-HIT, DIAMOND, MetaBAT2 |
| RNA-seq | FastQC, fastp, HISAT2, featureCounts, DESeq2, clusterProfiler |
| Statistics and figures | R (tidyverse, vegan, phyloseq, ggplot2, DESeq2) |
