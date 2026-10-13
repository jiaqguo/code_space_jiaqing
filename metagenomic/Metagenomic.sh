#!/bin/bash
#SBATCH -p bigmem
#SBATCH --job-name="coassemble_compare"
#SBATCH -t 168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mail-user=jiaqing.guo@okstate.edu
#SBATCH --mail-type=end

source /scratch/josh/software/miniconda3/bin/activate   #activate conda
conda activate metagenomic       #activate environment


#checkm2      :checkm2(v1.0.2)
#maxbin2      :maxbin2(v2.2.7)
#metabat2     :metabat2(v2.17)
#concoct      :concoct(v1.1.0) 
#metagenomic  :fastqc(v0.12.1), fastp(v0.23.4), multiqc(v1.22.3), bowtie2(v2.5.4),
#              samtools(v1.2.0), megahit(v1.2.9), prodigal(v2.6.3), cdhit(v4.8.1), diamond(v2.1.9),
#              fastqpair(v1.0) 





# quality check FastQC 

mkdir /scratch/josh/metagenomics-jing/results/fastqc_results/
indir=/scratch/josh/metagenomics-jing/rawdata/*fq.gz
outdir=/scratch/josh/metagenomics-jing/results/fastqc_results/
fastqc --threads 30 --outdir ${outdir} ${indir}
#conbine fastqc files  (v1.22.3)
mkdir /scratch/josh/metagenomics-jing/results/fastqc_results/multiqc/
outdir_mqc=/scratch/josh/metagenomics-jing/results/fastqc_results/multiqc/
multiqc ${outdir} --outdir ${outdir_mqc}


#read long: 50 bp 
# quality control   fastp   #2h
workdir=/scratch/josh/metagenomics-jing/rawdata/
mkdir /scratch/josh/metagenomics-jing/results/trim_data4/
mkdir /scratch/josh/metagenomics-jing/results/fastp_report4/
outdir=/scratch/josh/metagenomics-jing/results/trim_data4/
outdir_report=/scratch/josh/metagenomics-jing/results/fastp_report4/
cd ${workdir}
for i in *_1.fq.gz;  #define file name as variable
do
    i=${i/_1.fq.gz/}   
fastp --in1 ${i}_1.fq.gz \
      --in2 ${i}_2.fq.gz \
      --out1  ${outdir}${i}_1.clean.fq.gz \
      --out2  ${outdir}${i}_2.clean.fq.gz \
      --detect_adapter_for_pe \
      --length_required 50 \
      --cut_front \
      --cut_front_window_size 5 \
      --cut_front_mean_quality 30 \
      --cut_tail \
      --cut_tail_window_size 5 \
      --cut_tail_mean_quality 30 \
      --trim_poly_g \
      --average_qual 30 \
      --thread 1 \
      --correction \
      --html ${outdir_report}${i}_fastp.html
done

fastp --in1 NEs5_1.clean.fq.gz \
      --in2 NEs5_2.clean.fq.gz \
      --out1  NEs5_1.clean2.fq.gz \
      --out2  NEs5_2.clean2.fq.gz

# quality check again
#quality check FastQC (v0.12.1)
mkdir /scratch/josh/metagenomics-jing/results/fastqc_results_trim4/
indir=/scratch/josh/metagenomics-jing/results/trim_data4/*fq.gz
outdir=/scratch/josh/metagenomics-jing/results/fastqc_results_trim4/
fastqc --threads 30 --outdir ${outdir} ${indir}
#conbine fastqc files  (v1.22.3)
mkdir /scratch/josh/metagenomics-jing/results/fastqc_results_trim4/multiqc/
outdir_mqc=/scratch/josh/metagenomics-jing/results/fastqc_results_trim4/multiqc/
multiqc ${outdir} --outdir ${outdir_mqc}


# remove host genome bowtie2 (v2.5.4)  
#index genome using bowtie2  #2h
db_dir=/scratch/josh/metagenomic_hf_db/
host_fna=chicken.fna,maize.fna,soybean.fna
cd ${db_dir}
bowtie2-build ${host_fna} chick_food --threads 10

#remove host gene  #33.5h
mkdir /scratch/josh/metagenomics-jing/results/align_sam
mkdir /scratch/josh/metagenomics-jing/results/mvhost_data
db_index=/scratch/josh/metagenomic_hf_db/chick_food
data_dir=/scratch/josh/metagenomics-jing/results/trim_data4/
samout=/scratch/josh/metagenomics-jing/results/align_sam/
otudir=/scratch/josh/metagenomics-jing/results/mvhost_data/
cd ${data_dir}
for i in *_1.clean.fq.gz;
do 
    i=${i/_1.clean.fq.gz/}
    bowtie2 -x ${db_index} \
         -1 ${i}_1.clean.fq.gz \
         -2 ${i}_2.clean.fq.gz \
         -S ${samout}${i}.sam \
         --un-conc-gz ${otudir}${i} \
         --seed 07092024 \
         --threads 30 \
         --fr
done

grep "overall alignment rate" remove_host.txt | wc 

#modify file name
data_dir=/scratch/josh/metagenomics-jing/results/mvhost_data/
cd ${data_dir}
for i in *; 
do  
    if [[ "$i" == *".1" ]]; then
        cp ${i} ${i/.1/_1.fq.gz}
    elif [[ "$i" == *".2" ]]; then
        cp ${i} ${i/.2/_2.fq.gz}
    fi
done


# assemble
#Megahit (v1.2.9)
# k-min 21  k-max 141  k-step 12  [default]
# min-contig-len 200 [default]
# min-contig-len 500

#individual assemble  #30h
mkdir -p /scratch/josh/metagenomics-jing/results/assemble/megahit_individual/
data_dir=/scratch/josh/metagenomics-jing/results/mvhost_data/
outdir=/scratch/josh/metagenomics-jing/results/assemble/megahit_individual/
cd ${data_dir}
#loop run
for i in *_1.fq.gz;
do 
    i=${i/_1.fq.gz/}
    megahit \
        -1 ${i}_1.fq.gz \
        -2 ${i}_2.fq.gz \
        --k-min 21 \
        --k-max 141 \
        --k-step 12 \
        --min-contig-len 500 \
        --out-dir ${outdir}${i} \
        --out-prefix megahit_${i} \
        -t 30
done
#move individual assemble contigs to one folder
cd ${outdir}
for i in */;
do
    cd ${i}
    mv *.fa ../
    cd ..
done


#parallel run  #24h
mkdir -p /scratch/josh/metagenomics-jing/results/assemble/try/
data_dir=/scratch/josh/metagenomics-jing/results/mvhost_data/
cd ${data_dir}
ls *_1.fq.gz | parallel --jobs 2 -k '
    i={}
    i=${i/_1.fq.gz/}
    outdir=/scratch/josh/metagenomics-jing/results/assemble/try/
    megahit \
        -1 ${i}_1.fq.gz \
        -2 ${i}_2.fq.gz \
        --k-min 21 \
        --k-max 141 \
        --k-step 12 \
        --min-contig-len 500 \
        --out-dir ${outdir}${i} \
        --out-prefix megahit_${i} \
        -m 0.4 \
        -t 15
'

#co-assemble   #16.5h
indir=/scratch/josh/metagenomics-jing/results/mvhost_data/
outdir=/scratch/josh/metagenomics-jing/results/assemble/megahit_co/
cd ${indir}
R1=$(ls *_1* | tr "\n" "," | sed 's/,$//')
R2=$(ls *_2* | tr "\n" "," | sed 's/,$//')
megahit \
        -1 ${R1} \
        -2 ${R2} \
        --k-min 21 \
        --k-max 141 \
        --k-step 12 \
        --min-contig-len 500 \
        --out-dir ${outdir} \
        --out-prefix megahit_co \
        -t 30


# Gene pridiction
#Prodigal (v2.6.3)
# -d nucleotide sequences   -a protein translations  -s start/stop pairs along with score information

#co-assemble pridict    #80min
mkdir -p /scratch/josh/metagenomics-jing/results/prodigal/prodigal_co/
input=/scratch/josh/metagenomics-jing/results/assemble/megahit_co/megahit_co.contigs.fa
outdir=/scratch/josh/metagenomics-jing/results/prodigal/prodigal_co/
prodigal \
        -p meta \
        -i ${input} \
        -f gff \
        -o ${outdir}co.gff \
        -a ${outdir}co_protein_file.faa\
        -d ${outdir}co_mrna_file.fasta

#individual-assemble pridict
#loop run   #15.2h
mkdir -p /scratch/josh/metagenomics-jing/results/prodigal/prodigal_in/
indir=/scratch/josh/metagenomics-jing/results/assemble/megahit_individual/
outdir=/scratch/josh/metagenomics-jing/results/prodigal/prodigal_in/
cd ${indir}
for i in *fa;
do
    name=${i/megahit_/}
    name=${name/.contigs.fa/}
    prodigal \
        -p meta \
        -i ${i} \
        -f gff \
        -o ${outdir}in_${name}.gff \
        -a ${outdir}in_${name}_protein_file.faa\
        -d ${outdir}in_${name}_mrna_file.fasta
done

#parallel run  #49min
mkdir -p /scratch/josh/metagenomics-jing/results/prodigal/try2/
indir=/scratch/josh/metagenomics-jing/results/assemble/megahit_individual/
cd ${indir}
ls *fa | parallel --jobs 30 -k '
    i={}
    name=${i/megahit_/}
    name=${name/.contigs.fa/}
    outdir=/scratch/josh/metagenomics-jing/results/prodigal/try2/
    prodigal \
        -p meta \
        -i ${i} \
        -f gff \
        -o ${outdir}in_${name}.gff \
        -a ${outdir}in_${name}_protein_file.faa \
        -d ${outdir}in_${name}_mrna_file.fasta
'

# before de-redundancy
# combine all files together
# for mRNA
mkdir /scratch/josh/metagenomics-jing/results/cdhit/
codir=/scratch/josh/metagenomics-jing/results/prodigal/prodigal_co/
indir=/scratch/josh/metagenomics-jing/results/prodigal/try/
workdir=/scratch/josh/metagenomics-jing/results/cdhit/
workdir=/scratch/josh/metagenomics-jing/compar/
cp ${codir}*fasta ${workdir}
parallel --j 30 --bar cp {} distdir ::: ${workdir}/*

parallel --j 30 -k --bar cp -a ${indir}*fasta ::: ${workdir}
# for protein
cp ${codir}*faa ${workdir}
parallel --j 30 -k --bar cp -a ${indir}*faa ::: ${workdir}

#delete /n
#for mRNA
workdir=/scratch/josh/metagenomics-jing/results/cdhit/
cd ${workdir}
ls *_mrna_file.fasta | parallel --jobs 30 -k --bar '
    i={}
    name=${i/_mrna_file.fasta/}
    awk '"'"'
        BEGIN { RS=">"; FS="\n"; ORS="" } 
        NR > 1 { print ">" $1; print "\n";
            for (j=2; j<=NF; j++) print $j; print "\n" 
        }
    '"'"' ${i} > ${name}_remove.fasta
'
# for protein
ls *_protein_file.faa | parallel --jobs 30 -k --bar '
    i={}
    name=${i/_protein_file.faa/}
    awk '"'"'
        BEGIN { RS=">"; FS="\n"; ORS="" } 
        NR > 1 { print ">" $1; print "\n";
            for (j=2; j<=NF; j++) print $j; print "\n" 
        }
    '"'"' ${i} > ${name}_remove.faa
'

rm *rna* *protein* #delete orginal files

#filter mRNA less than 100bp
workdir=/scratch/josh/metagenomics-jing/results/cdhit/
cd "${workdir}"
ls *_remove.fasta | parallel --jobs 30 -k --bar '
    i={}
    name=${i/_remove.fasta/_filtered.fasta}
    awk '"'"'{
        if (NR % 2 == 1) {
            line1 = $0
        } else {
            if (length($0) >= 100) {
                print line1
                print $0
            }
        }
    }'"'"' ${i} > ${name}
'
#modify protein file according to mRNA file
ls *_filtered.fasta | parallel --jobs 30 -k --bar '
    gene_file={}
    protein_file=${gene_file/_filtered.fasta/_remove.faa}
    out_file=${gene_file/_filtered.fasta/_filtered.faa}
    awk '"'"'
        FNR==NR {
            if ($0 ~ /^>/) {
            names[$0] = 1
            }
            next
        }
        {if ($0 in names) {
            print $0
            getline
            print $0
            }
        }
    '"'"' ${gene_file} ${protein_file} > ${out_file}
'
rm *remove*  #delete unneeded

#add sample name to seqs names to avoid duplicate gene names
# for mRNA 
workdir=/scratch/josh/metagenomics-jing/results/cdhit/
cd ${workdir}
ls *_filtered.fasta | parallel --jobs 30 -k --bar '
    i={}
    name=${i/_filtered.fasta/}
    sed -i "s/>/>${name}_/" ${i}
'
# for protein
ls *_filtered.faa | parallel --jobs 30 -k --bar '
    i={}
    name=${i/_filtered.faa/}
    sed -i "s/>/>${name}_/" ${i}
'

#combine all files together
cat *.fasta > prodigal_all.fasta
cat *.faa > prodigal_all.faa
#rm *filter* 

# de-redundancy
#CD-Hit (v4.8.1)  #2.5h
#-c	sequence identity threshold [default 0.9]  -G 1:global sequence identity, 0:local sequence identity
#-M memory limit  -T number of threads  -n word_length
#-aS alignment coverage for the shorter sequence
#-g 1 more accurate
workdir=/scratch/josh/metagenomics-jing/results/cdhit/
workdir=/scratch/josh/metagenomics-jing/compar/
cd ${workdir}

cd-hit-est \
    -i prodigal_all.fasta \
    -o rep_all_mRNA.fasta \
    -c 0.95 \
    -g 1 \
    -G 0 \
    -M 0  \
    -T 0  \
    -aS 0.9


cd-hit-est \
    -i prodigal_individual.fasta \
    -o rep_individual_mRNA.fasta \
    -c 0.95 \
    -g 1 \
    -G 0 \
    -M 0  \
    -T 0  \
    -aS 0.9

cd-hit-est \
    -i prodigal_coassemble.fasta \
    -o rep_coassemble_mRNA.fasta \
    -c 0.95 \
    -g 1 \
    -G 0 \
    -M 0  \
    -T 0  \
    -aS 0.9

#get de-redundent protein
awk '
    FNR==NR {
        if ($0 ~ /^>/) {
        names[$0] = 1
        }
         next
    }
    {if ($0 in names) {
        print $0
        getline
        print $0
        }
    }
    ' rep_mRNA.fasta prodigal_all.faa > rep_protein.faa

#diamond (v2.1.9)





#binning
#filter contigs < 1500 bp
mkdir -p /scratch/josh/metagenomics-jing/results/bin/contigs/
#for co-assemle contigs
workdir=/scratch/josh/metagenomics-jing/results/assemble/megahit_co/
file=megahit_co.contigs.fa
outfile=/scratch/josh/metagenomics-jing/results/bin/contigs/co_filter.fa
cd ${workdir}
awk '{
        if (NR % 2 == 1) {
            line1 = $0
        } else {
            if (length($0) >= 1500) {
                print line1
                print $0
            }
        }
    }' ${file} > ${outfile}

#for individual assemble contigs
workdir=/scratch/josh/metagenomics-jing/results/assemble/megahit_individual/
cd ${workdir}
ls *.fa | parallel --jobs 30 -k --bar '
    i={}
    outdir=/scratch/josh/metagenomics-jing/results/bin/contigs/
    name=${i/megahit_/}
    name=${name/.contigs.fa/_filtered.fa}
    awk '"'"'{
        if (NR % 2 == 1) {
            line1 = $0
        } else {
            if (length($0) >= 1500) {
                print line1
                print $0
            }
        }
    }'"'"' ${i} > ${outdir}${name}
'

#get coverage table (map reads to co_assemble_filtered_contigs)
#index genome using bowtie2 
db=/scratch/josh/metagenomics-jing/results/bin/contigs/co_filter.fa
outdir=/scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/
cd ${outdir}
bowtie2-build ${db} co_filter --threads 30

#alignment
#parallel  #5.6h
mkdir -p /scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/samout
data_dir=/scratch/josh/metagenomics-jing/results/trim_data4/
cd ${data_dir}
ls *_1.clean.fq.gz | parallel --jobs 2 -k '
    db_index=/scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/co_filter
    samout=/scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/samout/
    R1={}
    R2=${R1/_1.clean.fq.gz/_2.clean.fq.gz}
    name=${R1/_1.clean.fq.gz/}
    bowtie2 -x ${db_index} \
         -1 ${R1} \
         -2 ${R2} \
         -S "${samout}${name}".sam \
         --seed 07212024 \
         --threads 15 \
         --fr
'

#loop #5.6h
mkdir -p /scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/samout2
data_dir=/scratch/josh/metagenomics-jing/results/trim_data4/
db_index=/scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/co_filter
samout=/scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/samout2/
cd ${data_dir}
for i in *_1.clean.fq.gz;
do 
    i=${i/_1.clean.fq.gz/}
    bowtie2 -x ${db_index} \
         -1 ${i}_1.clean.fq.gz \
         -2 ${i}_2.clean.fq.gz \
         -S ${samout}${i}.sam \
         --seed 07212024 \
         --threads 30 \
         --fr
done

#sam to bam  #74min
workdir=/scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/samout2/
cd ${workdir}
ls *sam | parallel --jobs 2 -k '
    i={}
    name=${i/.sam/.bam}
    samtools sort ${i} \
        -@ 15 \
        -o ${name}
'


#using metabat2
source /scratch/josh/software/miniconda3/bin/activate   #activate conda
conda activate metabat2      

#TNF (tetra-nucleotide frequency) binning  #30sample 6min
mkdir -p /scratch/josh/metagenomics-jing/results/bin/metabat2/bin_results/
datadir=/scratch/josh/metagenomics-jing/results/bin/contigs/
cd ${datadir}
ls *fa | parallel --jobs 2 -k '
    i={}
    name=${i/_filter.fa/}
    basename=metabat_${name}
    outdir=/scratch/josh/metagenomics-jing/results/bin/metabat2/bin_results/
    metabat2 \
        -i ${i} \
        -o ${outdir}${basename} \
        -t 15 \
        --seed  072124
'

#TNF (tetra-nucleotide frequency) and ABD (abundance) binning  #4min
outdir=/scratch/josh/metagenomics-jing/results/bin/metabat2/bin_results/
contigs_file=/scratch/josh/metagenomics-jing/results/bin/contigs/co_filter.fa
bam_dir=/scratch/josh/metagenomics-jing/results/bin/metabat2/coverage/samout2/
cd ${bam_dir}
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
metabat2 \
    -i ${contigs_file} \
    -a depth.txt \
    -o ${outdir}metabat_both \
    -t 30 \
    --seed  072124































