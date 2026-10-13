#qiime version: qiime2-amplicon-2023.9 

#download data

wd=/scratch/josh/repeat_guild

cd $wd
# download data
while IFS= read -r accession
do
    echo "Downloading $accession ...
    ----------------------------------"
    prefetch $accession
done < access_num.txt 

# file process
rm access_num.txt

for file in */;
do
    cd $file
    mv * ..
    cd ..
done

rm -r */

# uncompress data
fasterq-dump *

# delete unused files
rm *sra


# qiime2 part -------------------------------
#generate manifest file
wd="/scratch/josh/ziyu_test/data/"
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > temp.txt
for r1 in *_1*;
do
  id="${r1/_1.fastq/}"  
  r2="${r1/_1/_2}"
  echo -e "${id}\t${wd}${r1}\t${wd}${r2}" >> temp.txt
done
awk -F "\t" 'BEGIN {OFS="\t"} {gsub("_",".",$1); print}' temp.txt > manifest.txt && rm temp.txt



# create folder
results=/scratch/josh/ziyu_test/results/
scripts=/scratch/josh/ziyu_test/scripts/
classifier=/scratch/josh/workshop/classifier/
mkdir $results $scripts $classifier
cd $classifier


# Generate sbatch file and submit job.
script_dir=/scratch/josh/workshop/scripts/
cd $script_dir
vim qiime.sbatch
# 1. Press i
# 2. Copy all the code to the file
# 3. Type: :wq  to save and exit the file


# code for analysis----------------------------

#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name="qiime"
#SBATCH -t 20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end

source /scratch/josh/software/miniconda3/bin/activate   #activate conda
conda activate qiime2-amplicon-2025.4

results=/scratch/josh/workshop/results/
raw_dir=/scratch/josh/ziyu_test/rawdata/

cd $results

####import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${raw_dir}manifest.txt \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux.qza \
    --p-cores 10 \
    --p-overlap 3 \
    --p-front-f CCTAYGGGRBGCASCAG \
    --p-front-r GGACTACNNGGGTATCTAAT \
    --p-no-indels \
    --p-discard-untrimmed \
    --o-trimmed-sequences trimmed.qza 
#--p-overlap   [defult 3]
#--p-indels / --p-no-indels  [defult True]
#---p-discard-untrimmed / --p-no-discard-untrimmed [defult False]

#visualization 
qiime demux summarize \
    --i-data demux.qza \
    --o-visualization demux.qzv
qiime demux summarize \
    --i-data trimmed.qza \
    --o-visualization trimmed.qzv

####join the sequences
time qiime vsearch merge-pairs \
    --i-demultiplexed-seqs  trimmed.qza \
    --o-merged-sequences  joined.qza \
    --o-unmerged-sequences un_joined.qza \
    --p-minovlen 10 \
    --p-maxdiffs 10 \
    --p-threads 10
#--p-minovlen [defult 10]
#--p-maxdiffs [defult 10] 
#--p-threads [defult 10]

####qualify the sequences
time qiime quality-filter q-score \
    --i-demux joined.qza  \
    --p-min-quality 30 \
    --p-quality-window 5 \
    --p-min-length-fraction 1 \
    --o-filtered-sequences joined-filtered.qza \
    --o-filter-stats joined-filtered-stats.qza
#--p-quality-window [default 3]
#--p-min-length-fraction [default 0.75]

#visualization
qiime demux summarize \
    --i-data joined.qza \
    --o-visualization joined.qzv
qiime demux summarize \
    --i-data joined-filtered.qza \
    --o-visualization joined-filtered.qzv
qiime metadata tabulate \
    --m-input-file joined-filtered-stats.qza \
    --o-visualization  joined-filtered-stats.qzv

####denoise
time qiime deblur denoise-16S \
    --i-demultiplexed-seqs joined-filtered.qza \
    --p-trim-length 402 \
    --p-sample-stats  \
    --p-jobs-to-start 10 \
    --p-min-reads 10 \
    --p-min-size 1 \
    --o-table feature-table1.qza \
    --o-representative-sequences rep-seqs1.qza \
    --o-stats stats-deblur1.qza 
#--p-min-reads [defult 10]
#--p-min-size [defult 2]

#visualization
qiime feature-table tabulate-seqs \
    --i-data rep-seqs1.qza \
    --o-visualization rep-seqs1.qzv
qiime deblur visualize-stats \
  --i-deblur-stats stats-deblur1.qza \
  --o-visualization stats-deblur1.qzv
qiime feature-table summarize \
    --i-table feature-table1.qza \
    --o-visualization feature-table1.qzv 

###build the phylogenic tree
time qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza \
    --o-alignment seqs-aligned.qza \
    --o-masked-alignment rep-seqs_masked.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza \
    --p-n-threads auto




###taxonomic annotation
qiime feature-classifier classify-sklearn \
    --i-reads rep-seqs.qza \
    --i-classifier $classifier \
    --o-classification tax.qza

#data export
qiime tools export \
    --input-path feature-table.qza\
    --output-path ./
biom convert \
    -i feature-table.biom \
    -o feature_table.tsv \
    --to-tsv
qiime tools export \
    --input-path rep-seqs.qza \
    --output-path ./
qiime tools export \
    --input-path rooted-tree.qza \
    --output-path ./rooted_tree.nwk
qiime tools export \
    --input-path tax.qza \
    --output-path ./












#picrust2

#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name="picrust2"
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end
source /scratch/josh/software/miniconda3/bin/activate   #activate conda
conda activate picrust2

#after asv table pre-processing
picrust2_pipeline.py -s dna-sequences.fasta \
    -i feature_table_filter.tsv \
    -o picrust \
    -p 10 \
    --stratified



# rdp classifier

#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name="qiime"
#SBATCH -t 20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end
source /scratch/josh/software/miniconda3/bin/activate   #activate conda
conda activate qiime2-amplicon-2025.4

db=/scratch/josh/classifier/rdp_classifier.qza
result_dir=/scratch/josh/workshop/results
cd $result_dir


wget https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo19_QiimeFormat.zip
unzip RDPClassifier_16S_trainsetNo19_QiimeFormat.zip && rm RDPClassifier_16S_trainsetNo19_QiimeFormat.zip
mv RDPClassifier_16S_trainsetNo19_QiimeFormat/* ./ && rm -r RDPClassifier_16S_trainsetNo19_QiimeFormat

qiime tools import \
  --input-path RefOTUs.fa \
  --output-path rdp_dna.qza \
  --type 'FeatureData[Sequence]'

qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path  Ref_taxonomy.txt \
    --output-path  rdp_taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads rdp_dna.qza \
    --i-reference-taxonomy rdp_taxonomy.qza \
    --o-classifier rdp_classifier.qza





qiime feature-classifier classify-sklearn \
    --i-reads rep-seqs.qza \
    --i-classifier rdp_classifier.qza \
    --o-classification tax.qza

qiime tools export \
    --input-path tax.qza \
    --output-path ./








