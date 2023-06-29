## MAIN PIPELINE FILE

## This is the main wrapper file to run the EasySci-RNA computational processing pipeline. 
## It takes the raw FASTQ file folder, the sample ID file, the output folder, the STAR index folder, two GTF files (one with every genomic elements and one with only the exons), the file with the random hexamer and shortdT reverse transcription barcode matching, the number of cores and the flag indicating if single-end or paired-end processing will be used and produces the processed data under the "report/out/" folder.

## Inputs files and folders
fastq_folder=$1
sample_ID=$2
export output_folder=$3
index=$4
gtf_file=$5
gtf_file_exons=$6
RT_barcode_matching_file=$7
cores=$8
sequencing_type=$9

## Barcode information files
script_folder=$(dirname "$0")/script_folder/
## define the location of the ligation barcodes
ligation_barcode_file=$script_folder/barcode_files/ligation_barcodes.pickle2
## define the location of the RT barcodes
RT_barcode_file=$script_folder/barcode_files/RT_barcodes.pickle2
## define the location of the randomN RT barcodes
randomN_barcode_file=$script_folder/barcode_files/RT_randomN_barcodes.txt


## Change the name of the fastq files to standard format
echo "Changing the name of the fastq files..."
for sample in $(cat $sample_ID); do echo changing name $sample; mv $fastq_folder/*$sample*R1*.fastq.gz $fastq_folder/$sample.R1.fastq.gz; mv $fastq_folder/*$sample*R2*.fastq.gz $fastq_folder/$sample.R2.fastq.gz; mv $fastq_folder/*$sample*R3*.fastq.gz $fastq_folder/$sample.R3.fastq.gz; done
echo


#### Run the paired-end pipeline
if [ $sequencing_type == "paired-end" ]
then

## Barcode the reads
echo "Barcoding reads..."
script=$script_folder/barcoding_reads_paired.py
mkdir -p $output_folder/barcoded_fastqs
python $script $fastq_folder $sample_ID $output_folder/barcoded_fastqs $ligation_barcode_file $RT_barcode_file $cores $randomN_barcode_file
echo "Done barcoding reads"


## Trim the reads
echo "Start trimming..."
mkdir -p $output_folder/trimmed_fastqs/

trim() {
    echo "Trimming sample: $1"
    trim_galore --paired $output_folder/barcoded_fastqs/$1*R1*.gz $output_folder/barcoded_fastqs/$1*R2*.gz -a2 AAAAAAAA --stringency 3 -o $output_folder/trimmed_fastqs/
}

export -f trim

parallel -j $cores trim ::: $(cat $sample_ID)
echo "Done trimming..."
echo


## STAR alignment
echo "Start alignment using STAR..."
mkdir -p $output_folder/STAR_alignment
STAR --genomeDir $index --genomeLoad Remove

for sample in $(cat $sample_ID)
do 
echo "Aligning $sample..."
STAR --runThreadN $cores --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $output_folder/trimmed_fastqs/$sample*R1*gz $output_folder/trimmed_fastqs/$sample*R2*gz --outFileNamePrefix $output_folder/STAR_alignment/$sample --genomeLoad LoadAndKeep
done

STAR --genomeDir $index --genomeLoad Remove
echo "Done aligning"
echo


## Transfer STAR log files
mkdir -p $output_folder/report/Log_files
mv Log.out $output_folder/report/Log_files/
mv Log.progress.out $output_folder/report/Log_files/
mv Aligned.out.sam $output_folder/report/Log_files/
mv SJ.out.tab $output_folder/report/Log_files/
mv Log.final.out $output_folder/report/Log_files/


## Sorting and filtering SAM files
echo "Start sorting and filtering..."
mkdir -p $output_folder/filtered_sam

sort_filter() {
    echo "Sorting and filtering $1"
    samtools view -q 30 -f 2 -F 780 $output_folder/STAR_alignment/$1*.sam | sort -k1,1 -k3,3 -k4,4n > $output_folder/filtered_sam/$1.noheader.sam
    grep "@" $output_folder/STAR_alignment/$1*.sam > $output_folder/filtered_sam/$1.header.sam
    cat $output_folder/filtered_sam/$1.header.sam  $output_folder/filtered_sam/$1.noheader.sam > $output_folder/filtered_sam/$1.sam
    rm $output_folder/filtered_sam/$1.header.sam
    rm $output_folder/filtered_sam/$1.noheader.sam
    echo Filtering $1 done.
}
export -f sort_filter

parallel -j $cores sort_filter ::: $(cat $sample_ID)
echo "Done sorting and filtering"
echo


## Removing duplicates
echo "Start removing duplicates..."
mkdir -p $output_folder/duplicates_removed
script=$script_folder/duplicate_removal_paired.py
python $script $output_folder/filtered_sam/ $sample_ID $output_folder/duplicates_removed/ $cores
echo "Done removing duplicates."
echo


## Calculate read numbers along the pipeline
echo "Start calculating the reads number along the pipeline..."
mkdir -p $output_folder/report/read_num
echo sample,total reads,after filtering barcode,after trimming,uniquely aligned reads,After remove duplicates > $output_folder/report/read_num/read_number.csv
for sample in $(cat $sample_ID); do echo calculating $sample; echo $sample,$(expr $(zcat $fastq_folder/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $output_folder/barcoded_fastqs/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $output_folder/trimmed_fastqs/$sample*R2*.gz|wc -l) / 4),$(expr $(samtools view $output_folder/filtered_sam/$sample.sam|wc -l) / 2),$(expr $(samtools view $output_folder/duplicates_removed//$sample.sam|wc -l) / 2) >> $output_folder/report/read_num/read_number.csv; done
echo "Read number calculation is done."
echo


## Count the genes
echo "Start the gene count...."
mkdir -p $output_folder/report/Gene_count/ 
script=$script_folder/gene_counting_paired.py
python $script $gtf_file $output_folder/duplicates_removed/ $output_folder/report/Gene_count/ $sample_ID $cores $randomN_barcode_file
echo "Gene count done"
echo


## Post processing genes
echo "Start the gene post processing...."
script=$script_folder/post_processing_genes.py
mkdir -p $output_folder/report/out/genes/
python $script $output_folder/report/Gene_count/ $output_folder/report/out/genes/ $sample_ID $RT_barcode_matching_file
echo "Done with the gene post processing"
echo


## Count the exons
echo "Start the exon count...."
mkdir -p $output_folder/report/Exon_count/ 
script=$script_folder/exon_counting_paired.py
python $script $gtf_file_exons $output_folder/duplicates_removed/ $output_folder/report/Exon_count/ $sample_ID $cores
echo "Exon count done"
echo


## Post processing exons
echo "Start the exon post processing...."
script=$script_folder/post_processing_exons.py
mkdir -p $output_folder/report/out/exons/
python $script $output_folder/report/Exon_count/ $output_folder/report/out/exons/ $sample_ID $RT_barcode_matching_file $sequencing_type
echo "Done with the exon post processing"
echo


#### Run the single-end pipeline
else


## Barcode the reads
echo "Barcoding reads..."
script=$script_folder/barcoding_reads_single.py
mkdir -p $output_folder/barcoded_fastqs
python $script $fastq_folder $sample_ID $output_folder/barcoded_fastqs $ligation_barcode_file $RT_barcode_file $cores
echo "Done barcoding reads"


## Trim the reads
echo "Start trimming..."
mkdir -p $output_folder/trimmed_fastqs/

trim() {
    echo "Trimming sample: $1"
    trim_galore $output_folder/barcoded_fastqs/$1*R2*.gz -a AAAAAAAA --stringency 3 --three_prime_clip_R1 1 -o $output_folder/trimmed_fastqs/
}

export -f trim

parallel -j $cores trim ::: $(cat $sample_ID)
echo "Done trimming..."
echo


## STAR alignment
echo "Start alignment using STAR..."
mkdir -p $output_folder/STAR_alignment
STAR --genomeDir $index --genomeLoad Remove

for sample in $(cat $sample_ID)
do 
echo "Aligning $sample..."
STAR --runThreadN $cores --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $output_folder/trimmed_fastqs/$sample*R2*gz --outFileNamePrefix $output_folder/STAR_alignment/$sample --genomeLoad LoadAndKeep
done

STAR --genomeDir $index --genomeLoad Remove
echo "Done aligning"
echo


## Transfer STAR log files
mkdir -p $output_folder/report/Log_files
mv Log.out $output_folder/report/Log_files/
mv Log.progress.out $output_folder/report/Log_files/
mv Aligned.out.sam $output_folder/report/Log_files/
mv SJ.out.tab $output_folder/report/Log_files/
mv Log.final.out $output_folder/report/Log_files/


## Sorting and filtering SAM files
echo "Start sorting and filtering..."
mkdir -p $output_folder/filtered_sam

sort_filter() {
    echo "Sorting and filtering $1"
    samtools view -q 30 -F 772 $output_folder/STAR_alignment/$1*.sam | sort -k1,1 -k3,3 -k4,4n > $output_folder/filtered_sam/$1.noheader.sam
    grep "@" $output_folder/STAR_alignment/$1*.sam > $output_folder/filtered_sam/$1.header.sam
    cat $output_folder/filtered_sam/$1.header.sam  $output_folder/filtered_sam/$1.noheader.sam > $output_folder/filtered_sam/$1.sam
    rm $output_folder/filtered_sam/$1.header.sam
    rm $output_folder/filtered_sam/$1.noheader.sam
    echo Filtering $1 done.
}
export -f sort_filter

parallel -j $cores sort_filter ::: $(cat $sample_ID)
echo "Done sorting and filtering"
echo


## Removing duplicates
echo "Start removing duplicates..."
mkdir -p $output_folder/duplicates_removed
script=$script_folder/duplicate_removal_single.py
python $script $output_folder/filtered_sam/ $sample_ID $output_folder/duplicates_removed/ $cores
echo "Done removing duplicates."
echo


## Calculate read numbers along the pipeline
echo "Start calculating the reads number along the pipeline..."
mkdir -p $output_folder/report/read_num
echo sample,total reads,after filtering barcode,after trimming,uniquely aligned reads,After remove duplicates > $output_folder/report/read_num/read_number.csv
for sample in $(cat $sample_ID); do echo calculating $sample; echo $sample,$(expr $(zcat $fastq_folder/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $output_folder/barcoded_fastqs/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $output_folder/trimmed_fastqs/$sample*R2*.gz|wc -l) / 4),$(samtools view $output_folder/filtered_sam/$sample.sam|wc -l),$(samtools view $output_folder/duplicates_removed//$sample.sam|wc -l) >> $output_folder/report/read_num/read_number.csv; done
echo "Read number calculation is done."
echo


## Count the genes
echo "Start the gene count...."
mkdir -p $output_folder/report/Gene_count/ 
script=$script_folder/gene_counting_single.py
python $script $gtf_file $output_folder/duplicates_removed/ $output_folder/report/Gene_count/ $sample_ID $cores $randomN_barcode_file
echo "Gene count done"
echo


## Post processing genes
echo "Start the gene post processing...."
script=$script_folder/post_processing_genes.py
mkdir -p $output_folder/report/out/genes/
python $script $output_folder/report/Gene_count/ $output_folder/report/out/genes/ $sample_ID $RT_barcode_matching_file
echo "Done with the gene post processing"
echo


## Count the exons
echo "Start the exon count...."
mkdir -p $output_folder/report/Exon_count/ 
script=$script_folder/exon_counting_single.py
python $script $gtf_file_exons $output_folder/duplicates_removed/ $output_folder/report/Exon_count/ $sample_ID $cores
echo "Exon count done"
echo


## Post processing exons
echo "Start the exon post processing...."
script=$script_folder/post_processing_exons.py
mkdir -p $output_folder/report/out/exons/
python $script $output_folder/report/Exon_count/ $output_folder/report/out/exons/ $sample_ID $RT_barcode_matching_file $sequencing_type
echo "Done with the exon post processing"
echo


fi


## Remove intermediate folders
rm -r $output_folder/barcoded_fastqs/
rm -r $output_folder/trimmed_fastqs/
rm -r $output_folder/STAR_alignment/
rm -r $output_folder/filtered_sam/

echo "EasySci-RNA pipeline completed"


