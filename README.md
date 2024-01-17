# EasySci
## Computational pipeline to process EasySci-RNA data.

This pipeline is for processing EasySci-RNA datasets and takes the raw FASTQ files as inputs and outputs the cell/gene (or exon) matrix. The pipeline consists of the following steps: barcode extraction and matching; trimming of the adaptor and polyA sequences, alignment to the reference genome; filtering of the low-quality alignments; PCR duplicate removal; gene and exon expression counting per single cell.

Please note: This computational pipeline is an updated version of the pipeline used to the process the EasySci dataset in our publication. The updates include making the pipeline faster, saving space on disk by eliminating single-cell SAM file generation and an integrated post-processing step that merges the reads from the two different reverse transcription primer layers. Small differences are expected compared to the original pipeline, particularly the updated version is more sensitive and will recover slightly more cells. To reproduce the results in our publication, pleaser refer to the original pipeline we used to process the data: https://zenodo.org/records/8395492.

## Dependencies:

star: version 2.7.10b  
samtools: version 1.17  
trim-galore: version 0.6.10   
python: version 3.10.11  
htseq: version 2.0.3  
pandas: version 1.5.2   
numpy: version 1.23.5   
scipy: version 1.10.0  
biopython: version 1.81   
parallel: version 20230522  



## Pre-processing: to generate the input FASTQ files, follow this demultiplexing tutorial:

Demultiplexing for EasySci-RNA uses only the P7 barcodes, because in EasySci-RNA the P5 barcode is added by ligation and therefore every possible P5-P7 combination is present in the final data. This would result in a very high number of files if both barcodes would be used for demultiplexing. The demultiplexing is done with Illumina’s bcl2fastq software with the following settings:

bcl2fastq --runfolder-dir INPUT_FOLDER(sequencer generated files) -o OUTPUT_FOLDER --sample-sheet SAMPLE_SHEET--reports-dir OUTPUT_FOLDER/report --barcode-mismatches 1 --create-fastq-for-index-reads --no-lane-splitting --use-bases-mask Y*,I*,Y*,Y* --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0

This demultiplexing generates 4 files:  
•	R1: Read 1  
•	R2: P5 barcode  
•	R3: Read 2  
•	I1: P7 barcode  


## Running the EasySci-RNA computational pipeline:

Run the following code to use the EasySci-RNA computational pipeline after the demultiplexing step:

bash path_to_script/EasySci_pipeline.sh $fastq_folder $sample_ID $output_folder $STAR_index $gtf_file $gtf_file_exons $RT_barcode_matching_file $cores $sequencing_type

The following input parameters need to be set and must follow the order presented below:

•	fastq_folder: Folder with the input FASTQ files.

•	sample_ID: Sample ID file, name of the demultiplexed files, one sample name per row without the R1/R1/R3 or fastq ending.

•	output_folder: Output folder for the intermediary and final files.

•	STAR_index: STAR index folder for alignment.

•	gtf_file: GTF file for the gene counting step.

•	gtf_file_exons: GTF file for the exon counting step, this file can be created from the above GTF file by subsetting for only the exonic regions, removing the redundant exons, changing the feature column from ‘exon’ to ‘gene’ and formatting the final column. An example of this file for both human (GRCh38) and mouse (GRCm39) can be found under ‘/script_folder/reference_files/’.

•	RT_barcode_matching_file: This file is used to match the random hexamer and shortdT reverse transcription reads from the same well. The following file can be used if barcodes from the reference publication were used, and reverse transcription wells received primers from the corresponding primer wells. ‘/script_folder/barcode_files/RT_barcode_matching.txt’

•	cores: Number of cores used during the computational pipeline.

•	sequencing_type: ‘paired-end’ or ‘single-end’, this decision should be made based on the length of Read1. If ~60 bp were sequenced from Read1, then the ‘single-end’ flag should be used, while if ~100 bp were sequenced from Read1, then use the ‘paired-end’ flag. Generally, to make the decision, the length of Read1 after removing the barcodes and polyA tail should be considered (18 + 15 = 33 bp). If the sequence length after removing 33 bp are sufficient to reliable align to the genome, then the paired-end pipeline can be used, which utilizes both Read 1 and Read2 for gene and exon counting. 


## Output of the EasySci-RNA pipeline:

The final output files can be found here:

‘output_folder/report/out/genes/’  
‘output_folder/report/out/exons/’  

•	cell_annotation.csv: cell annotation csv file  
•	gene_annotation.csv/exon_annotation.csv: gene or exon annotation csv file  
•	expression_matrix.mtx: gene or exon x cell expression matrix in sparse matrix format (matrix market format)  

The cells are not filtered by any criteria, but during the computational pipeline the reads originating from the shortdT and random hexamer RT primers of the same cells are merged. The pipeline assumes that the barcoded oligos from the reference publication were used in the reverse transcription and ligation steps.


## Reference:

Andras Sziraki, Ziyu Lu, Jasper Lee, Gabor Banyai, Sonya Anderson, Abdulraouf Abdulraouf, Eli Metzner, Andrew Liao, Jason Banfelder, Alexander Epstein, Chloe Schaefer, Zihan Xu, Zehao Zhang, Li Gan, Peter T. Nelson, Wei Zhou, Junyue Cao. A global view of aging and Alzheimer’s pathogenesis-associated cell population dynamics and molecular signatures in human and mouse brains. Nature Genetics 55, 2104–2116 (2023). https://doi.org/10.1038/s41588-023-01572-y

## Credits

The EasySci-RNA computational pipeline was developed by Andras Sziraki, M.D. at the Cao Lab, The Rockefeller University.
