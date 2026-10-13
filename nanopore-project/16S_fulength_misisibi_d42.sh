#!/bin/bash
#SBATCH -p bullet
#SBATCH --job-name="dorado"
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:rtx6000:2
#SBATCH --mem=100G
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end


#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name="file modify"
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end


########## software and version ##########
dorado v0.9.1
samtools v1.21
cutadapt v5.0
cdhit v4.8.1
minimap2 v2.28
medaka v2.0.1
vsearch v2.30.0
emu v3.5.1
########## software and version ##########

source /scratch/josh/software/miniconda3/bin/activate 
conda activate pod5

# dorado
model=/scratch/josh/software/dorado-0.9.1-linux-x64/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
rawdata=/scratch/josh/16_fulength/pod5


cd /scratch/josh/16_fulength
dorado basecaller $model $rawdata \
    --min-qscore 20 \
    --kit-name SQK-NBD114-96 > simplex22.bam
#  --barcode-both-ends 
#  --emit-moves  polish need (assemble)

# check process
dorado summary simplex6.bam > summary_simplex.tsv
# print unclassified percentage
num1=$(grep unclassified summary_simplex6.tsv | wc -l)
num6=$(cat summary_simplex6.tsv | wc -l)
result=$(awk "BEGIN {print $num1 / $num6}")
echo "Result: $result"


# demux after basecaller (simplex)
dorado demux --output-dir demux6 --no-classify simplex6.bam

# merger two runs reads
samtools merge barcode01.bam *barcode01.bam -@ $(nproc)

# check barcodes
samtools view barcode01.bam | grep  AAGAAAGTTGTCGGTGTCTTTGTG | wc -l


# bam to fastq file
for file in ls *bam;
  do name=${file/.bam/}
  samtools fastq $file > ${name}.fastq
done

ls *bam | parallel --jobs 20 -k --bar '
    file={}
    name=${file/.bam/}
  samtools fastq $file > ${name}.fastq
'


# downstream analysis

##################################################################
####################    Cut adapter     ##########################
##################################################################
conda activate nanopore
# cutadapter for 16S-primer removing and length filtering
fwd_primer=CGGTTACCTTGTTACGACTT
rev_primer=AGAGTTTGATCCTGGCTCAG
RC_fwd_primer=AAGTCGTAACAAGGTAACCG
RC_rev_primer=CTGAGCCAGGATCAAACTCT

cd /scratch/josh/16_fulength/rawfq/
outdir=/scratch/josh/16_fulength/results/cutprimer/
for file in *.fastq; do
  name=${file/.fastq/}
  cutadapt \
      ${file} \
      -g ^${fwd_primer} \
      -a ${RC_rev_primer}$ \
      -n 2 \
      -j $(nproc) \
      -o ${outdir}cutadapt_tmp_${name}.fastq \
   && cutadapt \
      ${outdir}cutadapt_tmp_${name}.fastq \
      -g ^${rev_primer} \
      -a ${RC_fwd_primer}$ \
      -n 2 \
      -j $(nproc)\
      --minimum-length 1300 \
      --maximum-length 1950 \
      -o ${outdir}cutadapt_${name}.fastq
done
cd /scratch/josh/16_fulength/results/cutprimer/
rm *tmp*.fastq


# remove non-16S primer sequences （use 16Sreads_filter.sh script)
outdir="/scratch/josh/16_fulength/results/filtered_fq"
input_base="/scratch/josh/16_fulength/results/cutprimer/cutadapt_"
refer_base="/scratch/josh/16_fulength/rawfq/"
cd /scratch/josh/16_fulength/rawfq/
ls *fastq > /scratch/josh/16_fulength/ids.txt
# Read the IDs from the file
while read -r id; do
    # Set the input and reference file paths for each ID
    input_fq="${input_base}${id}"
    refer_fq="${refer_base}${id}"

    # Run the 16S filtering script with the dynamic filenames
    sh 16Sreads_filter.sh --input "$input_fq" --reference "$refer_fq" --outdir "$outdir"
done < ids.txt


# get reads stats 
beforedir=/scratch/josh/16_fulength/rawfq/
afterdir=/scratch/josh/16_fulength/results/cutprimer/
filterdir=/scratch/josh/16_fulength/results/filtered_fq/
outdir=/scratch/josh/16_fulength/results/
seqkit stats -a ${beforedir}*fastq >> raw_stats.txt
seqkit stats -a ${afterdir}*fastq >> cutadapt_stats.txt
seqkit stats -a ${filterdir}*fastq >> cut_filter_stats.txt





##################################################################
###############   EMU Taxonomy Classification  ###################
##################################################################
conda activate nanopore

# download database
pip install osfclient
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
tar -xvf emu.tar


#!/bin/bash
#SBATCH -p bigmem
#SBATCH --job-name="emu_build"
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end

source /scratch/josh/software/miniconda3/bin/activate 
conda activate nanopore

# build self-database
# python

#------------------------
#!/usr/bin/env python3
import gzip
import sys

if len(sys.argv) != 4:
    print("Usage: python build_seq2tax.py database.fasta nucl_gb.accession2taxid.gz seq2tax.map.tsv")
    sys.exit(1)

fasta_file = sys.argv[1]
acc2taxid_file = sys.argv[2]
output_file = sys.argv[3]

# Step 1: read header in FASTA
accessions = []
with open(fasta_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            # get accession 
            acc = line[1:].strip().split()[0]
            accessions.append(acc)

accessions = set(accessions)  # dereplication

# Step 2: read NCBI accession2taxid.gz create dictionary
acc2tax = {}
with gzip.open(acc2taxid_file, "rt") as f:
    next(f)  # skip header
    for line in f:
        cols = line.strip().split("\t")
        acc = cols[1]
        taxid = cols[2]
        if acc in accessions:
            acc2tax[acc] = taxid

# Step 3: output seq2tax.map.tsv
with open(output_file, "w") as out:
    for acc in accessions:
        if acc in acc2tax:
            out.write(f"{acc}\t{acc2tax[acc]}\n")
        else:
            sys.stderr.write(f"WARNING: {acc} not found in accession2taxid\n")


python3 build_seq2tax.py bacteria.16SrRNA.fna nucl_gb.accession2taxid.gz seq2tax.map.tsv
#------------------------

awk 'BEGIN{OFS="\t"; print "id","lineage"} 
{
    split($0, arr, "\t")
    id = arr[1]
    n = split(arr[2], tax, ";")
    lineage = ""
    for (i=n; i>=1; i--) {
        gsub(/^ +| +$/, "", tax[i]) 
        if (lineage=="") {
            lineage = tax[i]
        } else {
            lineage = lineage ";" tax[i]
        }
    }
    print id, lineage
}' 666.txt > taxonomy.tsv

awk 'NR==FNR{map[$1]=$2; next} 
    FNR==1{print $0; next} 
    ($1 in map){$1=map[$1]} 
    {print $0}' seq2tax.map.tsv taxonomy.tsv > updated.tax.tsv


# -*- coding: utf-8 -*-
import csv

# Input files
map_file = "merged.dmp"           # First file (old_id -> new_id)
seq_file = "seq2tax.map.tsv"      # Second file
out_file = "seq2tax.map.replaced.tsv"  # Output file

# 1. Read mapping (old_id -> new_id)
id_map = {}
with open(map_file, "r", encoding="utf-8") as f:
     for line in f:
        # Remove leading/trailing whitespace
        line = line.strip()
        # Split by '|', strip whitespace from each part
        parts = [p.strip() for p in line.split("|")]
        if len(parts) >= 2:
            old_id = parts[0]
            new_id = parts[1]
            id_map[old_id] = new_id

print(f"Loaded {len(id_map)} mapping entries")

# 2. Process the second file
with open(seq_file, "r", encoding="utf-8") as fin, open(out_file, "w", encoding="utf-8") as fout:
    for line in fin:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            acc = parts[0]
            taxid = parts[1]
            if taxid in id_map:  # Replace taxid if it exists in the mapping
                taxid = id_map[taxid]
            fout.write(f"{acc}\t{taxid}\n")

print(f"Processing finished. Result saved to {out_file}")





#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name="emu_build"
#SBATCH -t 30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end

source /scratch/josh/software/miniconda3/bin/activate 
conda activate nanopore

cd /scratch/josh/16_fulength/results/db_ncbi/
emu build-database real_emu_db \
    --sequences bacteria.16SrRNA.fna \
    --seq2tax seq2tax.map.replaced.tsv \
    --ncbi-taxonomy /scratch/josh/16_fulength/results/db_ncbi/



# run

db=/scratch/josh/16_fulength/results/db_ncbi/real_emu_db/
ref_seq=/scratch/josh/16_fulength/results/filtered_fq/
outdir=/scratch/josh/16_fulength/results/emu
for file in ${ref_seq}*.fastq; do
    emu abundance $file \
    --type map-ont \
    --db $db \
    --output-dir $outdir \
    --threads $(nproc)
done



# epi2me results


#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name="epi2me"
#SBATCH -t 20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end

module load apptainer

export NXF_HOME=/scratch/josh/nextflow_home
export APPTAINER_CACHEDIR=/scratch/josh/singularity_cache
export APPTAINER_TMPDIR=/scratch/josh/tmpdir
export TMPDIR=/scratch/josh/tmpdir
export NXF_SINGULARITY_CACHEDIR=/scratch/josh/singularity_cache
export NXF_SINGULARITY_CLI=apptainer

fq_dir=/scratch/josh/16_fulength/results/fq_epi2me/
nextflow run /scratch/josh/nextflow_home/assets/epi2me-labs/wf-16s \
    -profile singularity \
	--fastq $fq_dir \
    --sample_sheet /scratch/josh/16_fulength/results/barcodes.csv \
	--taxonomic_rank S \
    --keep_bam \
    --minimap2_by_reference \
    --abundance_threshold 0 \
    --include_read_assignments \
    --output_unclassified \
    --threads 20

# for unclassified reads
fq_dir=/scratch/josh/16_fulength/results/fq_unclassified
nextflow run /scratch/josh/nextflow_home/assets/epi2me-labs/wf-16s \
    -profile singularity \
	--fastq $fq_dir \
    --sample_sheet /scratch/josh/16_fulength/results/samples.csv \
    --classifier kraken2 \
	--taxonomic_rank S \
    --abundance_threshold 0 \
    --include_read_assignments \
    --output_unclassified \
    --threads 10


nextflow run /scratch/josh/nextflow_home/assets/epi2me-labs/wf-16s \
    -profile singularity \
	--fastq $fq_dir \
    --classifier kraken2 \
	--taxonomic_rank S \
    --abundance_threshold 0 \
    --include_read_assignments \
    --output_unclassified \
    --threads 10



for file in $(ls); do 
    name=${file/.unclassified.fq/}
    mkdir $name
    mv $file $name
done


echo "barcode,alias" > samples.csv
for s in sample*; do
    echo "$s,$s" >> samples.csv
done


for file in *; do
  name=${file/sample/barcode}
  mv $file $name
done

sed -i -E 's/^sample([0-9]+)/barcode\1/' samples.csv
# ()捕获组  \1 第一个捕获组内容

for file in *; do
  cd $file
  cp * ../
  cd ../
done




#teset
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-16s/wf-16s-demo.tar.gz
tar -xzvf wf-16s-demo.tar.gz

nextflow run /scratch/josh/nextflow_home/assets/epi2me-labs/wf-16s \
	--fastq /scratch/josh/16_fulength/results/test/wf-16s-demo/test_data \
	-profile singularity




#######################      not use       #######################

##################################################################
################      cluster process        #####################
##################################################################
conda activate nanopore    # use CD-HIT

# -c	sequence identity threshold [default 0.9]  
# -n word_length [default 10]
# -g 1 more accurate [default 0]
# -d length of description in .clstr file [default 20]
# -M memory limit  
# -T number of threads 

file=/scratch/josh/16_fulength/all.fa
outfile=/scratch/josh/16_fulength/results/cluster/cluster.fa
cd-hit-est \
    -i "$file" \
    -o "$outfile" \
    -c 0.97 \
    -n 11 \
    -g 1 \
    -d 0 \
    -M 0  \
    -T 0 

##################################################################
################       polish process        #####################
##################################################################

# Generate clust and reads information -----------

# python script
import os
import re

# Set file paths
cluster_info_file = "/scratch/josh/16_fulength/results/cluster/cluster.fa.clstr"
cluster_seq_file = "/scratch/josh/16_fulength/results/cluster/cluster.fa"
reads_file = "/scratch/josh/16_fulength/all.fa"
output_dir = "/scratch/josh/16_fulength/results/cluster/cluster_seq_files/"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Generate cluster_map
cluster_map = {}

with open(cluster_info_file, 'r') as f:
    current_cluster = ""
    for line in f:
        if line.startswith(">Cluster"):
            current_cluster = line.strip().replace(" ", "").replace(">", "")
        elif "..." in line:
            match = re.search(r">(.*?)[.]", line)
            read_id = ">" + match.group(1)
            cluster_map[read_id] = current_cluster

print("[INFO] mapping file done")

# Generate cluster_seq files
def process_sequences(input_file, output_suffix):
    with open(input_file, 'r') as f:
        current_cluster = None
        next_line = False
        for line in f:
            if line.startswith(">"):
                read_id = line.strip()
                if read_id in cluster_map:
                    current_cluster = cluster_map[read_id]
                    output_file = os.path.join(output_dir, f"{current_cluster}{output_suffix}.fa")
                    with open(output_file, 'a') as out_f:
                        out_f.write(f"{read_id}\n")
                    next_line = True
                else:
                    next_line = False
            elif next_line:
                with open(output_file, 'a') as out_f:
                    out_f.write(line)

# Process cluster_seq_file
process_sequences(cluster_seq_file, "")

print("[INFO] cluster sequence done")

# Process reads_file
process_sequences(reads_file, "_reads")

print("[INFO] reads sequence done")
# python script done-----------------


# Medaka for consensus reads (optional for dorado)
conda activate nanopore

mkdir -p /scratch/josh/16_fulength/results/cluster/bam_files
mkdir -p /scratch/josh/16_fulength/results/cluster/polish_hdf
mkdir -p /scratch/josh/16_fulength/results/cluster/polish_results
cd /scratch/josh/16_fulength/results/cluster/cluster_seq_files

parallel -j 12 -a ../job_list.txt \
    'mini_align -i Cluster{1}_reads.fa -r Cluster{1}.fa -m \
        -p ../bam_files/aligned_reads{1} -t 1'
echo "[INFO] alignment done"

parallel --halt now,fail=1 -j 12 -a ../job_list.txt \
    'medaka inference ../bam_files/aligned_reads{1}.bam ../polish_hdf/output{1}.hdf \
        --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0:consensus'
echo "[INFO] inference done"

parallel --halt now,fail=1 -j 12 -a ../job_list.txt \
    'medaka sequence ../polish_hdf/output{1}.hdf ../polish_results/polished{1}.fa'
echo "[INFO] all done"

# check model to use (optional)
medaka tools resolve_model --auto_model consensus <input.bam/input.fastq>

##################################################################
################            denoise          #####################
##################################################################

# add abundance annotations-------------
# python scripts
from tqdm import tqdm

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values.
    
    :param file_path: Path to the FASTA file.
    :return: Dictionary of sequences.
    """
    with open(file_path, 'r') as file:
        sequences = {}
        current_id = None
        for line in tqdm(file, desc="Parsing FASTA file"):
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]  # Remove '>' and take the sequence ID
                sequences[current_id] = ''
            else:
                sequences[current_id] += line  # Append sequence data
    return sequences

def parse_clusters(file_path):
    """
    Parses a cluster file and returns a dictionary with cluster numbers as keys and lists of sequence IDs as values.
    
    :param file_path: Path to the cluster file.
    :return: Dictionary of clusters.
    """
    clusters = {}
    current_cluster = None
    with open(file_path, 'r') as file:
        for line in tqdm(file, desc="Parsing Cluster file"):
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = line  # Extract cluster
                clusters[current_cluster] = []
            else:
                seq_id = line.split('>')[1].split('...')[0].strip()  # Extract sequence ID
                clusters[current_cluster].append(seq_id)
    return clusters

def process_sequences(sequences, clusters):
    """
    Processes sequences by updating their IDs to include cluster information.
    
    :param sequences: Dictionary of sequences.
    :param clusters: Dictionary of clusters.
    :return: Dictionary of processed sequences.
    """
    # Create a lookup dictionary for sequence IDs to their cluster and cluster sizes
    id_to_cluster = {}
    cluster_sizes = {cluster: len(clusters[cluster]) for cluster in clusters}
    for cluster, ids in clusters.items():
        for seq_id in ids:
            id_to_cluster[seq_id] = cluster
    #
    processed_sequences = {}
    for seq_id, sequence in tqdm(sequences.items(), desc="Processing sequences"):
        if seq_id in id_to_cluster:
            cluster = id_to_cluster[seq_id]
            new_id = f">{seq_id};size={cluster_sizes[cluster]};otu={seq_id}"
            processed_sequences[new_id] = sequence
    return processed_sequences

def write_fasta(sequences, output_path):
    """
    Writes sequences to a FASTA file.
    
    :param sequences: Dictionary of sequences.
    :param output_path: Path to the output FASTA file.
    """
    with open(output_path, 'w') as file:
        for seq_id, sequence in tqdm(sequences.items(), desc="Writing to FASTA file"):
            file.write(f"{seq_id}\n")
            file.write(f"{sequence}\n")

def main(fasta_path, clusters_path, output_path):
    """
    Main function to process FASTA and cluster files and produce an updated FASTA file.
    
    :param fasta_path: Path to the input FASTA file.
    :param clusters_path: Path to the cluster file.
    :param output_path: Path to the output FASTA file.
    """
    #
    sequences = parse_fasta(fasta_path)
    clusters = parse_clusters(clusters_path)
    if sequences is None or clusters is None:
            print("\033[91mError: Failed to parse input files\033[0m")
            return
    
    processed_sequences = process_sequences(sequences, clusters)
    if processed_sequences is None:
        print("\033[91mError: 'processed_sequences' is None\033[0m")
        return
    write_fasta(processed_sequences, output_path)

if __name__ == "__main__":
    from datetime import datetime, timedelta
    start_time = datetime.now()
    
    fasta_path = 'cluster.fa'  # Path to the input FASTA file
    clusters_path = 'cluster.fa.clstr'  # Path to the cluster file
    output_path = 'cluster_size.fa'  # Path to the output FASTA file
    
    print(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    main(fasta_path, clusters_path, output_path)
    
    end_time = datetime.now()
    duration = end_time - start_time
    print(f"\nScript finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total duration: {str(timedelta(seconds=duration.seconds))}")
# python scripts done----------------------

# Sort and remove singletons
conda activate vsearch
echo Sort and remove singletons
echo

vsearch \
    --sortbysize cluster_size.fa \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --minsize 2 \
    --output sorted.fasta

echo
echo Non-singleton clusters: $(grep -c "^>" sorted.fasta)

# add sample annotations
# python scripts-----------------------
import os
import glob

def add_filename_to_header(input_fastq):
    filename = os.path.basename(input_fastq).split('_')[0]
    output_fastq = f"{filename}_name.fastq"
    
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        for line in infile:
            if line.startswith('@'):
                outfile.write(f"{line.strip()};sample={filename}\n")
            else:
                outfile.write(line)

fastq_files = glob.glob('*fastq')            
for fastq_file in fastq_files:
    add_filename_to_header(fastq_file)

# fastq to fa
import os

def fastq_to_fasta(input_fastq):
    base_name = os.path.splitext(input_fastq)[0]
    output_fasta = f"{base_name}.fa"
    #
    with open(input_fastq, 'r') as infile, open(output_fasta, 'w') as outfile:
        while True:
            header = infile.readline().strip()
            if not header:
                break
            sequence = infile.readline().strip()
            infile.readline()  # Skip the '+' line
            infile.readline()  # Skip the quality score line
            #
            fasta_header = f">{header[1:]}"
            outfile.write(f"{fasta_header}\n{sequence}\n")

if __name__ == "__main__":
    input_fastq = "all_name.fastq"
    fastq_to_fasta(input_fastq)
# python scripts done----------------------

# chimera removal
source /scratch/josh/software/miniconda3/bin/activate 
conda activate vsearch

echo 
echo De novo chimera detection
echo
vsearch --uchime3_denovo sorted.fasta \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --nonchimeras denovo.nonchimeras.fasta \

echo
echo Unique sequences after de novo chimera detection: $(grep -c "^>" denovo.nonchimeras.fasta)
echo
echo Reference chimera detection
echo

# refer_db.fa: 16S sequences from NCBI
vsearch --uchime_ref denovo.nonchimeras.fasta \
        --threads $(nproc) \
        --db /scratch/josh/16_fulength/refer_db.fa \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --dbmask none \
        --nonchimeras nonchimeras.fasta

echo
echo Unique sequences after reference-based chimera detection: $(grep -c "^>" nonchimeras.fasta)
echo
echo
echo Map sequences to OTUs by searching
echo

# relabel OTU
# some clusters were removed, therefore map again
input=/scratch/josh/16_fulength/results/cluster97/nonchimeras.fa
ref_db=/scratch/josh/16_fulength/results/cluster97/otus.fasta
reads=/scratch/josh/16_fulength/all_name.fa
output=/scratch/josh/16_fulength/results/cluster97/otutab.txt

vsearch --fastx_filter $input \
        --threads $(nproc) \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel OTU_ \
        --fastaout $ref_db

vsearch --usearch_global $reads \
        --threads $(nproc) \
        --db $ref_db \
        --id 0.97 \
        --strand both \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --dbmask none \
        --otutabout $output

echo OTU table done


# taxonomy analysis 
wd=/scratch/josh/classifier
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads $wd/rdp_dna.qza \
    --i-reference-taxonomy $wd/rdp_taxonomy.qza \
    --o-classifier $wd/rdp_classifier.qza

# import reads
dir=/scratch/josh/16_fulength/results/cluster97
qiime tools import \
  --input-path  ${dir}/otus.fasta \
  --output-path  ${dir}/dna.qza \
  --type 'FeatureData[Sequence]'

# taxonomy annotation kmer
dir=/scratch/josh/16_fulength/results/cluster97
classifier=/scratch/josh/classifier/rdp_classifier.qza
qiime feature-classifier classify-sklearn \
    --i-reads ${dir}/dna.qza \
    --i-classifier ${classifier} \
    --p-confidence 0.97 \
    --o-classification ${dir}/taxonomy_kmer.qza

qiime tools export \
    --input-path ${dir}/taxonomy_kmer.qza \
    --output-path ${dir}/taxonomy_kmer


# 2 database
# RDP: taxonomy_vsearch.qza blast6.qza taxonomy_blast.qza  blast_blast6.qza 
# NCBI: ncbi_taxonomy_vsearch.qza ncbi_blast6.qza ncbi_taxonomy_blast.qza ncbi_blast_blast6.qza

# taxonomy annotation vsearch
db_dir=/scratch/josh/ncbi_db
wk_dir=/scratch/josh/16_fulength/results/cluster97

qiime feature-classifier classify-consensus-vsearch \
    --i-query $wk_dir/dna.qza \
    --i-reference-reads $db_dir/ncbi_dna.qza \
    --i-reference-taxonomy $db_dir/ncbi_taxonomy.qza \
    --p-perc-identity 0.97 \
    --p-query-cov 0.99 \
    --p-threads 10 \
    --o-classification $wk_dir/ncbi_taxonomy_vsearch.qza \
    --o-search-results  $wk_dir/ncbi_blast6.qza

qiime tools export \
    --input-path $wk_dir/ncbi_taxonomy_vsearch.qza \
    --output-path $wk_dir/ncbi_taxonomy_vsearch

qiime tools export \
    --input-path  $wk_dir/ncbi_blast6.qza \
    --output-path  $wk_dir/ncbi_blast6



# taxonomy annotation blast
db_dir=/scratch/josh/ncbi_db
wk_dir=/scratch/josh/16_fulength/results/cluster97


qiime feature-classifier classify-consensus-vsearch \
    --i-query $wk_dir/dna.qza \
    --i-reference-reads $db_dir/ncbi_dna.qza \
    --i-reference-taxonomy $db_dir/ncbi_taxonomy.qza \
    --p-perc-identity 0.97 \
    --p-query-cov 0.99 \
    --p-threads 10 \
    --o-classification $wk_dir/ncbi_taxonomy_blast.qza \
    --o-search-results  $wk_dir/ncbi_blast_blast6.qza

qiime tools export \
    --input-path $wk_dir/ncbi_taxonomy_blast.qza \
    --output-path $wk_dir/ncbi_taxonomy_blast

qiime tools export \
    --input-path  $wk_dir/ncbi_blast_blast6.qza \
    --output-path  $wk_dir/ncbi_blast_blast6








#---------------------------------------------------------#
#                 process information                     #
#---------------------------------------------------------#

# TAAGGTTAA  before: adapters
# SQK-NBD114-96  barcode
# CACCCCAGTCATCTGCCCT  start of reads
# fwdprimer CGGTTACCTTGTTACGACTT
# revprimer AGAGTTTGATCCTGGCTCAG

# basecall and demux information
# single end barcode no-classify
945341 unclassifed
11669966 total   8.46%

# after cut adapter, filter and cluster, 99% cluster, -0 not accurate 6h
9390158 reads
618654  clusters
Non-singleton clusters: 75706
Remove chimera clusters: 54466
Reads retained:

# after cut adapter, filter and cluster, 97% cluster, -1 accurate 2d2h40min
9390158 reads  
127300  clusters
Non-singleton clusters: 48079
Remove chimera clusters: 33823
Reads retained:


#---------------------------------------------------------#
#                MAY used for duplex                      #
#---------------------------------------------------------#

# basecall your reads
# demultiplex reads
# extract read IDs from demultiplexed reads
# cleanup read ID files
# use pod5 tools and read ID files to demultiplex pod5 files by barcode
# split pod5 file
# duplex

# split reads by channel
rawdata=/scratch/josh/16_fulength/pod5
pod5 view $rawdata --include "read_id, channel" --output summary.tsv
less -S summary.tsv
pod5 subset $rawdata --summary summary.tsv --columns channel --output split_pod

model=/scratch/josh/software/dorado-0.9.1-linux-x64/dna_r10.4.1_e8.6_400bps_sup@v5.0.0
rawdata=/scratch/josh/16_fulength/pod5

dorado basecaller sup -r pod5_directory/  --kit-name SQK-NBD114-64 > basecalls.bam

dorado demux -t $(nproc) --output-dir basecalls_separated --no-classify basecalls.bam

for file in basecalls_separated/*.bam;
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  samtools view -@ $(nproc) -h $file | cut -f 1 > $id.readids.txt; 
  done

for file in *readids.txt; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  grep -v "@" $file > $id.clean.txt; 
  done

mkdir pod5_by_barcode

for file in *clean.txt; 
  do 
  id=$(echo "$file" | grep -oP 'barcode\d+'); 
  pod5 filter -r pod5_directory/ -i $file -t $(nproc) --missing-ok --duplicate-ok --output pod5_by_barcode/$id.pod5; 
  done


#---------------------------------------------------------#
#                  PRONAME pipeline                       #
#---------------------------------------------------------#
# set up image
# change cache directory
mkdir -p /scratch/josh/.apptainer
export APPTAINER_CACHEDIR=/scratch/josh/.apptainer  # run every time  or ln link
# load module
module load apptainer
# pull image
export APPTAINER_DOCKER_USERNAME=jiaqingguo
export APPTAINER_DOCKER_PASSWORD=GJQwasd1634
apptainer pull /scratch/josh/software/proname1.1.0.sif docker://benn888/proname:v1.1.0-amd64
# check image
apptainer inspect --environment proname1.1.0.sif


# run command image
export APPTAINER_CACHEDIR=/scratch/josh/.apptainer 
module load apptainer
image=/scratch/josh/software/proname1.1.0.sif
myenvs=/scratch/josh/software/myenvs
rawfq=/scratch/josh/16_fulength/rawfq

# remove host gene?
apptainer exec --env-file ${myenvs} --bind ${rawfq} $image proname_import \
  --inputpath ${rawfq} \
  --threads $(nproc) \
  --duplex no \
  --trimadapters no \
  --trimprimers yes \
  --fwdprimer CGGTTACCTTGTTACGACTT \
  --revprimer AGAGTTTGATCCTGGCTCAG \
  --noscatterplot yes

# check primer
grep AGAGTTTGATCCTGGCTCAG reads_simplex.fastq | wc -l

# length filter
apptainer exec --env-file ${myenvs} $image proname_filter \
  --datatype simplex \
  --filtminlen 1300 \
  --filtmaxlen 1950 \
  --filtminqual 0 \
  --threads $(nproc) \
  --inputpath ${rawfq} \
  --noscatterplot yes


# cluster reads
apptainer exec --env-file ${myenvs} $image proname_refine \
  --clusterid 0.99 \
  --clusterthreads $(nproc) \
  --medakathreads $(nproc) \
  --inputpath ${rawfq} \
  --medakamodel r1041_e86_400bps_sup_v5.0.0 \
  --chimeradb /opt/db/rEGEN-B/rEGEN-B_sequences.fasta \
  --medakabatchsize 150 \
  --qiime2import no


apptainer exec --env-file ${myenvs} $image proname_taxonomy


# run apptainer image interatively
export APPTAINER_CACHEDIR=/scratch/josh/.apptainer 
module load apptainer
mount_file=/scratch/josh/16_fulength/results
image=/scratch/josh/software/qiime2024.10.sif
myenvs=/scratch/josh/software/myenvs
apptainer shell --bind ${mount_file} --env-file ${myenvs} ${image}
apptainer run --bind ${mount_file} ${image}
PATH=/opt/conda/envs/qiime2-amplicon-2024.5/bin:$PATH






#---------------------------------------------------------#
#                        not used                         #
#---------------------------------------------------------#
# minimap2 remove host genome：  chicken,maize,soybean
host_genome1=/scratch/josh/chick_genome/GCA_000002315.5_GRCg6a_genomic.fna
host_genome2=/scratch/josh/chick_genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
host_genome3=/scratch/josh/chick_genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna
minimap2 -t $(nproc) -ax map-ont ${host_genome1} cutadapt_output.fastq | \
    samtools view -@ $(nproc) -b -f 4 | samtools fastq -@ $(nproc) > temp1.fastq
minimap2 -t $(nproc) -ax map-ont ${host_genome2} temp1.fastq | \
    samtools view -@ $(nproc) -b -f 4 | samtools fastq -@ $(nproc) > temp2.fastq
minimap2 -t $(nproc) -ax map-ont ${host_genome3} temp2.fastq | \
    samtools view -@ $(nproc) -b -f 4 | samtools fastq -@ $(nproc) > non_host.fastq && rm temp*fastq

seqkit stats -a non_host.fastq


# dorado polish
# Align the reads using dorado aligner, sort and index
source /scratch/josh/software/miniconda3/bin/activate 
conda activate samtools

cd /scratch/josh/16_fulength/results/cluster/cluster_seq_files/
seq 0 618653 > ../job_list.txt
parallel -j 4 -a ../job_list.txt \
    'dorado aligner Cluster{1}.fa Cluster{1}_reads.fa | \
     samtools sort > ../bam_files/aligned_reads{1}.bam && 
     samtools index ../bam_files/aligned_reads{1}.bam'

echo "[INFO] alignment done"

# add basecall model info as RG header in bam files
@HD     VN:1.6  SO:unknown
@RG     DS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0

# Call consensus  (GPU faster)
parallel -j 4 -a ../job_list.txt \
    'dorado polish ../bam_files/aligned_reads{1}.bam Cluster{1}.fa > ../polish_cluster/polished_Cluster{1}.fa'









##################################################################
###########          some check code              ################
##################################################################
# extract size information
# python
import csv
import re
# Function to extract size values from FASTA headers and write them to a CSV file
def extract_sizes_from_fasta(fasta_file, output_csv):
    """
    Extracts size values from the headers of a FASTA file and writes them to a CSV file.
    
    :param fasta_file: Path to the input FASTA file.
    :param output_csv: Path to the output CSV file.
    """
    # Define a regex pattern to match the size value in the FASTA header
    size_pattern = re.compile(r";size=(\d+)")
    #
    # List to store extracted size values
    sizes = []
    #
    # Read the FASTA file and process each line
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith(">"):  # Process only header lines
                match = size_pattern.search(line)
                if match:
                    sizes.append(match.group(1))  # Extract and store the size value
    #
    # Write the extracted size values to a CSV file
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["size"])  # CSV file header
        for size in sizes:
            writer.writerow([size])  # Write each size value to the CSV file

# Call the function with example data and specify the output CSV file name
extract_sizes_from_fasta('polish_cluster_size_sorted.fa', 'sizes.csv')


/scratch/josh/16_fulength/results/filtered_fq

# Generate reads file before polish
def parse_fastq(file_path):
    with open(file_path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            yield header, sequence, plus, quality

def find_common_heads(file1_path, file2_path, output_path):
    file2_headers = set()
    for header, _, _, _ in parse_fastq(file2_path):
        file2_headers.add(header)
    #
    with open(output_path, 'w') as output_file:
        for header, sequence, plus, quality in parse_fastq(file1_path):
            if header in file2_headers:
                output_file.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
    #  

if __name__ == "__main__":
    file1_path = "/scratch/josh/16_fulength/all_seq.fastq"
    file2_path = "/scratch/josh/16_fulength/results/cluster97/polish_cluster.fastq"
    output_path = "/scratch/josh/16_fulength/results/cluster97/before_polish.fastq"
    find_common_heads(file1_path, file2_path, output_path)

# filter quality
chopper --threads $(nproc) --quality 30 --input all_seq.fastq > all_seq_30.fastq

# check quality
NanoPlot -t $(nproc) --fastq before_polish.fastq -o polish_before
NanoPlot -t $(nproc) --fastq polish_cluster.fastq -o polish_after






