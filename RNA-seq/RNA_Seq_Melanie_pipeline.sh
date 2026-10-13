
#!/bin/bash
#SBATCH -p bigmem
#SBATCH -t 168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end

source /scratch/josh/software/miniconda3/bin/activate   #activate conda
conda activate rna-seq        #activate rna-seq environment

#quality check FastQC (v0.12.1)
fastqc --threads 30 --outdir /scratch/josh/Melanie/DCA_RNA_Seq/result/fastqc /scratch/josh/Melanie/DCA_RNA_Seq/rawdata/*fq.gz
#conbine fastqc files  (v1.22.3)
multiqc /scratch/josh/Melanie/DCA_RNA_Seq/result/fastqc


#read long: 50 bp for small RNA
#quality control   fastp (v0.23.4)
for i in $(ls /scratch/josh/Melanie/DCA_RNA_Seq/rawdata/*_1.fq.gz)  #define file name as variable
do
    i=${i/_1.fq.gz/}   #${x1/x2/x3} replace x2 as x3 in x1
    i=${i/\/scratch\/josh\/Melanie\/DCA_RNA_Seq\/rawdata\//}
fastp --in1 /scratch/josh/Melanie/DCA_RNA_Seq/rawdata/${i}_1.fq.gz \
      --in2 /scratch/josh/Melanie/DCA_RNA_Seq/rawdata/${i}_2.fq.gz \
      --out1  /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/${i}_1.clean.fq.gz \
      --out2  /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/${i}_2.clean.fq.gz \
      --length_required 50 \
      --adapter_sequence AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
      --adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
      --qualified_quality_phred 30 \
      --trim_front1=12 \
      --cut_right \
      --cut_window_size 5 \
      --cut_mean_quality 30 \
      --thread 10 \
      --correction \
      --html /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/fastp.html \
      --json /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/fastp.json
done

#--------------------------------------------------------------------------------
fastp --in1 /scratch/josh/Melanie/DCA_RNA_Seq/rawdata/*1.fq.gz \
      --in2 /scratch/josh/Melanie/DCA_RNA_Seq/rawdata/*2.fq.gz \
      --out1  /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/ \
      --out2  /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/ \
      --dont_overwrite \     #do not cover files
      --length_required 50 \    #min length
      --adapter_sequence = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
      --adapter_sequence_r2 = GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
      --qualified_quality_phred 30 \  #base quality
      --trim_front1=12 \     #cut read1 front 12 bases
      --trim_front2=12 \     #cut read2 front 12 bases
      --detect_adapter_for_pe \   #enable adapter sequence auto-detection （For PE data, the adapter sequence auto-detection is disabled by default since the adapters can be trimmed by overlap analysis.）
      --cut_right \   #cut from front to tail, drop base and the right in the window
      --cut_window_size 5 \   #window size
      --cut_mean_quality 30 \  #mean quality of the bases in the window
      --thread 10 \
      --correction \   #enable the base correction by overlap analysis
      --json \  #output json file
      --html   #output html
#--------------------------------------------------------------------------------

#quality check after trim
mkdir /scratch/josh/Melanie/DCA_RNA_Seq/result/fastqc_trim
fastqc --threads 10 --outdir /scratch/josh/Melanie/DCA_RNA_Seq/result/fastqc_trim /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/*fq.gz
#conbine fastqc files
multiqc --outdir /scratch/josh/Melanie/DCA_RNA_Seq/result/fastqc_trim \
                 /scratch/josh/Melanie/DCA_RNA_Seq/result/fastqc_trim



##read allignment

#make exon, splicesite file   Hisat2 (v2.2.1) 
hisat2_extract_splice_sites.py /scratch/josh/chicken_refer_genome/genome.gtf > /scratch/josh/chicken_refer_genome/genome.ss
hisat2_extract_exons.py /scratch/josh/chicken_refer_genome/genome.gtf > /scratch/josh/chicken_refer_genome/genome.exon

#build HFM index
hisat2-build -p 30 \
             --seed 021224 \
             --exon /scratch/josh/chicken_refer_genome/genome.exon \
             --ss /scratch/josh/chicken_refer_genome/genome.ss \
             -f /scratch/josh/chicken_refer_genome/genome.fna \
             /scratch/josh/chicken_refer_genome/chicken

#read allign to index
mkdir /scratch/josh/Melanie/DCA_RNA_Seq/result/allign
for i in $(ls /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/*_1.clean.fq.gz)  #define file name as variable
do
    i=${i/_1.clean.fq.gz/}   #${x1/x2/x3} replace x2 as x3 in x1
    i=${i/\/scratch\/josh\/Melanie\/DCA_RNA_Seq\/result\/trim\//}
    hisat2 -x /scratch/josh/chicken_refer_genome/chicken \
           -1 /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/${i}_1.clean.fq.gz \
           -2 /scratch/josh/Melanie/DCA_RNA_Seq/result/trim/${i}_2.clean.fq.gz \
           -S /scratch/josh/Melanie/DCA_RNA_Seq/result/allign/allign_${i}.sam \
           --seed 022024 \
           --phred33 \
           -k 2 \
           --threads 30
done


#check the progress
grep -o "overall alignment rate" slurm-3201875.out | wc
# -o show only the part of a line matching PATTERN

#---------------------------------------------
hisat2 -x /scratch/josh/chicken_refer_genome/chicken \  #basename for indexed genome
       -1 \  #pair-end one
       -2 \  ##pair-end two
       -S \  #Output alignment to file (SAM format) 
       --phred33 \
       --seed 022024 \
       --phred33 \
       -k 5 \  #search for 5 multimapping
       --dta \  #output alignments suitable for transcriptome assembly
       --rna-strandness RF \   # Specify read strand-specific information, forward or reverse
       --threads
#---------------------------------------------


#---------------------------------------
#read quantification, just show, don't need
flattenGTF GTF.featureType = "exon" \
           GTF.attrType = "gene_id" \
           method = "merge" \
           -a genome.gtf \
           -o genome.saf
#check gene number          
awk -F'\t' '{
     {print $1;}
}' genome.saf | uniq | wc
#---------------------------------------

#read summarization   featureCounts (v2.0.6) 
featureCounts -T 30 \
              -p \
              --countReadPairs \
              -C \
              -s 0 \
              -a /scratch/josh/chicken_refer_genome/genome.gtf \
              -o /scratch/josh/Melanie/DCA_RNA_Seq/result/counts.txt \
              /scratch/josh/Melanie/DCA_RNA_Seq/result/allign/*.sam


#------------------------------------------------------------
featureCounts -T 30 \  #thread
              -p \    #(isPairedEnd) Specify that input data contain paired-end reads.  
              --countReadPairs \
              #-M \   #(countMultiMappingReads) If specified, multi-mapping reads/fragments will be counted.
              #--fraction \ #each alignment carries 1/x count
              #-B  \  #(requireBothEndsMapped) If specified, only fragments that have both ends successfully aligned will be considered for summarization.
              -C \  #(countChimericFragments) If specified, the chimeric fragments (those fragments that have their two ends aligned to different chromosomes) will NOT be counted.
              -t exon \  #(GTF.featureType) Specify the feature type. [default: exon]
              -g gene_id \    #Specify the attribute type used to group features (eg. exons) into meta-features (eg. genes)
              #-f  \   #If specified, read summarization will be performed at featurelevel [default:meta-feature level]
              -s 0 \  #(isStrandSpecific) [default: 0] 0:unstranded; 1:stranded; 2:reversely stranded
              -F \  #specify annotation file type (GTF or SAF, [default: GTF])
              -a /scratch/josh/chicken_refer_genome/genome.gff \
              -o /scratch/josh/Melanie/DCA_RNA_Seq/result/counts.txt \
              /scratch/josh/Melanie/DCA_RNA_Seq/result/allign.*.sam
#----------------------------------------------------------------------------







