#!/bin/bash
##### help #######
usage() {
        cat <<-EOF
Usage: $0 [-w working_directory] [-d database_directory] [-c classifier] [-h]
Options:
    -w   Work directory, where is the rep-seqs (modify your file name as rep-seqs.qza)
    -d   Directory of the database
    -c   Specific name of your RDP classifier. (exp. /scratch/josh/classifier/rdp_classifier.qza)
    -h   Display help information
EOF
} 

while getopts "w:d:c:h" opt; do
    case "${opt}" in
        w)
            wd=${OPTARG}
            ;;
        d)
            database_dir=${OPTARG}
            ;;
        c)
            classifier=${OPTARG}
            ;;
        h)
            usage
            exit 0
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

#check if all the options exist
if [ -z "${wd}" ] || [ -z "${database_dir}" ] || [ -z "${classifier}" ]; then
    echo -e "[Error]: Missing required arguments."
    usage
    exit 1
fi

#add_slash
if [ "${wd: -1}" != "/" ]; then
        wd="${wd}/"
fi

if [ "${database_dir: -1}" != "/" ]; then
        database_dir="${database_dir}/"
fi


##### main script #####
#annotation
cd $wd
qiime feature-classifier classify-consensus-vsearch \
   --i-query  ${wd}rep-seqs.qza \
   --i-reference-reads  ${database_dir}self_dna.qza \
   --i-reference-taxonomy  ${database_dir}self_taxonomy.qza \
   --p-maxaccepts 'all' \
   --p-top-hits-only \
   --p-perc-identity 1 \
   --p-query-cov 1 \
   --p-min-consensus 1 \
   --p-strand 'both' \
   --p-threads 10 \
   --o-classification ${wd}classification.qza \
   --o-search-results ${wd}hits.qza

if [ $? -ne 0 ]; then
    exit 1
fi

qiime tools export \
    --input-path  ${wd}classification.qza \
    --output-path  ${wd}classification
qiime tools export \
    --input-path  ${wd}hits.qza \
    --output-path  ${wd}hits


###### Taxonomy check begin ####
#Run the following code if you need to look for annotation < 7 levels feature and print feature ID
#If you need to look for annotation < 7 levels feature, then find and run related code on the bottom of the document first, then go ahead to run the next code.
row_num_tax=$(cat classification/taxonomy.tsv | wc -l)
row_num_hit=$(($(cat hits/blast6.tsv | wc -l) + 1))
if [ ${row_num_tax} -eq ${row_num_hit} ]; then
    echo -e "[Info] Same row number" # 34m indicate blue
else
    echo -e "[WARNING] Different row number" # 33m indicate yellow
fi
###### Taxonomy check over ####


###Looking for EzBioCloud assigned features
grep -v "Unassigned" classification/taxonomy.tsv | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,"EzBioCloud"}' | sed '1d' > ezbcloud.txt

###Looking for Unassigned feature and print feature ID
grep "Unassigned" classification/taxonomy.tsv | awk -F '\t' '{print $1}' > temp.txt
###Looking for the fasta seqs of unassigned seqs
grep -A 1 -f temp.txt dna-sequences.fasta | grep -v "^--$" > unassigned.fasta && rm temp.txt
# grep -A 1 -f names.txt sequences.fasta：find seqs' name from names.txt in sequences.fasta, and also print the next row. 
# grep -v "^--$": -v means print rows not include ^--$. Last stemp will output these ^--$ between two seqs

###Annotation for Unassigned feature using RDP database
qiime tools import \
  --input-path unassigned.fasta \
  --output-path unassigned.qza \
  --type 'FeatureData[Sequence]'

qiime feature-classifier classify-sklearn \
    --i-reads unassigned.qza \
    --i-classifier ${classifier} \
    --p-confidence 0.97 \
    --o-classification unassigned-tax.qza
qiime tools export \
    --input-path unassigned-tax.qza \
    --output-path unassigned-tax


#add annotation source
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,"RDP"}' unassigned-tax/taxonomy.tsv | sed '1d' | awk -F'\t' '{
    split($2, fields, ";")
    output = ""
    for (i=1; i<=7; i++) {
        if (fields[i] != "") {
            output = output (i == 1 ? "" : ";") fields[i]
        } else {
            output = output (i == 1 ? "" : ";") fields[i-1]
            fields[i] = fields[i-1]
        }
    }
    print $1 "\t" output "\t" $3
}'  > rdp.txt

###Combine the taxonomy from EzBioCloud and RDP
cat ezbcloud.txt rdp.txt > taxonomy.txt
#move unused file
mkdir -p unused_file
mv hits.qza classification.qza hits classification unassigned* ezbcloud.txt rdp.txt ./unused_file
