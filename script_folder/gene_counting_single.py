'''
This script takes the files from the duplicate removed folder, extracts the single cells and then counts the gene expressions per cell for the single-end pipeline.
'''

import collections
import numpy as np
import pandas as pd
from multiprocessing import Pool
import multiprocessing
import HTSeq
from functools import partial
import sys

## Function to read in GTF file and to send reads to be counted paralell
def EasySci_count_parallel(gtf_file, input_folder, output_folder, sample_ID_file, core_number, randomN_barcode_file):
    # read in the gtf file, and then construct the genome interval for exons, genes, and gene end dictionary
    gtf_file = HTSeq.GFF_Reader(gtf_file, end_included=True)
    gene_annotat_file = output_folder + "/gene_name_annotate.csv"  
    gene_annotat = open(gene_annotat_file, "w")
    
    exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    genes = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    gene_end = {}
    exon_n = 0
    gene_n = 0
    transcript_n = 0
    gene_count = 0
    print("Start generating exon genomic arrays....")
    print("Start generating gene genomic arrays....")
    print("Start calculating transcript end of genes....")

    for feature in gtf_file:
        if feature.type == "exon":
            exon_n += 1
            exons[ feature.iv ] += feature.attr["gene_id"]
        elif feature.type == "gene":
            gene_n +=1
            genes[ feature.iv ] += feature.attr["gene_id"]
            gene_count += 1
            
            message = (feature.attr["gene_id"] + "," + feature.attr["gene_type"] + "," 
                       + "exon" + "," + feature.attr["gene_name"] + "," + str(gene_count) + "\n")
            
            gene_annotat.write(message)
            
            gene_count += 1
            
            message = (feature.attr["gene_id"] + "_intron" + "," + feature.attr["gene_type"] + "," 
                       + "intron" + "," + feature.attr["gene_name"] + "_intron" + "," + str(gene_count) + "\n")
            gene_annotat.write(message)
            
        elif feature.type == "transcript":
            transcript_n += 1
            if feature.attr["gene_id"] in gene_end.keys():
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)
            else:
                gene_end[ feature.attr["gene_id"] ] = set()
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)

    print("Detected gene number: ", gene_n)
    print("Detected transcript number: ", transcript_n)
    print("Detected exon number: ", exon_n)
    
    gene_annotat.close()
    
    gene_annotat = pd.read_csv(gene_annotat_file, header=None)
    gene_annotat.index =  gene_annotat[0]
    
    sample_ID = list(pd.read_csv(sample_ID_file, header=None)[0])
        
    # generate the randomN barcode list:
    barcodes = open(randomN_barcode_file, "rb")
    with barcodes as f:
        randomN_barcodes = f.read().splitlines()
    barcodes.close()

    # Parallel processing
    p = Pool(processes = int(core_number))
    func = partial(EasySci_count, input_folder, output_folder , exons, genes, gene_end, gene_annotat, randomN_barcodes)
    result = p.map(func, sample_ID)
    p.close()
    p.join()

    print("All analysis done~")

    
## Function to separate, count and save reads from different cells
def EasySci_count(input_folder, output_folder, exons, genes, gene_end, gene_annotat, randomN_barcodes, sample):
    input_file_name = input_folder + "/" + sample + ".sam"
    almnt_file = HTSeq.SAM_Reader(input_file_name)
    pre_barcode = 'EMPTY'
    cell_ID = 1
    counts = collections.Counter()
    count_output_name = output_folder + "/" + sample + ".count"
    count_output = open(count_output_name, "w")
    cell_annotat_file = output_folder + "/" + sample + "_cell_annotate.csv"
    cell_output = open(cell_annotat_file, "w")
  
    for alnmt in almnt_file:
        
        barcode = alnmt.read.name.split(',')[0]

        if barcode == pre_barcode:
            counts = count_read(alnmt, genes, exons, gene_end, randomN_barcodes, barcode, counts)
        else:
            if len(counts) == 0:
                counts = count_read(alnmt, genes, exons, gene_end, randomN_barcodes, barcode, counts)
                pre_barcode = barcode
            else:
                unmatched_reads = 0
                ambiguous_reads = 0
                all_reads = 0
                for gene in counts:
                    if (gene in ["_unmapped", "_ambiguous", "_ambiguous_intron", "_no_feature"]):
                        if gene in ["_no_feature"]:
                            unmatched_reads += counts[gene]
                        else:
                            ambiguous_reads += counts[gene]
                    else:
                        line = str(gene_annotat.loc[gene,4]) + "," + str(cell_ID) + "," + str(counts[gene]) + "\n"
                        count_output.write(line)
                        all_reads += counts[gene]

                unmatched_rate = unmatched_reads / (unmatched_reads + all_reads + ambiguous_reads)
                line_cell = sample + "." + pre_barcode + "," + str(unmatched_rate) + "\n"
                cell_output.write(line_cell)

                counts = collections.Counter()
                cell_ID += 1         
                counts = count_read(alnmt, genes, exons, gene_end, randomN_barcodes, barcode, counts)
                pre_barcode = barcode

    unmatched_reads = 0
    ambiguous_reads = 0
    all_reads = 0
    for gene in counts:
        if (gene in ["_unmapped", "_ambiguous", "_ambiguous_intron", "_no_feature"]):
            if gene in ["_no_feature"]:
                unmatched_reads += counts[gene]
            else:
                ambiguous_reads += counts[gene]
        else:
            line = str(gene_annotat.loc[gene,4]) + "," + str(cell_ID) + "," + str(counts[gene]) + "\n"
            count_output.write(line)
            all_reads += counts[gene]

    unmatched_rate = unmatched_reads / (unmatched_reads + all_reads + ambiguous_reads)
    line_cell = sample + "." + pre_barcode + "," + str(unmatched_rate) + "\n"
    cell_output.write(line_cell)

    count_output.close()
    cell_output.close()
    
    return 0     

    
## Function to assign one read pair to a gene
def count_read(alnmt, genes, exons, gene_end, randomN_barcodes, barcode, counts):
 
    sample_barcode = barcode[10:20]    
    RT_barcode = "shortDT"
    if sample_barcode in randomN_barcodes:
        RT_barcode = "randomN"
    
    
    if not alnmt.aligned:
        counts["_unmapped"] += 1      
    elif alnmt.iv.chrom not in genes.chrom_vectors:
        counts["_unmapped"] += 1
    else:  
        # First check the intersectin with exons
        gene_id_intersect = set()
        gene_id_combine = set()
        inter_count = 0
        for cigop in alnmt.cigar:
            if cigop.type != "M":
                continue

            for iv,val in exons[cigop.ref_iv].steps():
                gene_id_combine |= val
                if inter_count == 0:
                    gene_id_intersect |= val
                    inter_count += 1
                else:
                    gene_id_intersect &= val

        # first check the intersection set
        if len(gene_id_intersect) == 1:
            gene_id = list(gene_id_intersect)[0]
            counts[gene_id] += 1
        elif len(gene_id_intersect) > 1:
            gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end, RT_barcode)
            counts[gene_id] += 1
        else:
            # if there no intersection match, then find the union sets
            if len(gene_id_combine) == 1:
                gene_id = list(gene_id_combine)[0]
                counts[gene_id] += 1
            elif len(gene_id_combine) > 1:
                gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end, RT_barcode)
                counts[gene_id] += 1
            else:
                # if there is no intersection match or union match, then search for genes to find the intronic match
                gene_id_intersect = set()
                gene_id_combine = set()
                inter_count = 0
                for cigop in alnmt.cigar:
                    if cigop.type != "M":
                        continue
                    for iv,val in genes[cigop.ref_iv].steps():
                        gene_id_combine |= val
                        if inter_count == 0:
                            gene_id_intersect |= val
                            inter_count += 1
                        else:
                            gene_id_intersect &= val
                                                    
                ## If there is no intersection match or union match after the intronic search, test if there is an overlap in the 3' end flanking regions or on the opposite strand.    
                exonintron = "_intron"
                if len(gene_id_intersect) == 0 and len(gene_id_combine) == 0:
                    gene_id_combine, gene_id_intersect = test_flanking_region(alnmt, genes)
                if len(gene_id_intersect) == 0 and len(gene_id_combine) == 0:
                    gene_id_combine, gene_id_intersect = test_opposite_strand(alnmt, exons, 1)
                    if len(gene_id_intersect) > 0 or len(gene_id_combine) > 0:
                        exonintron = ""
                if len(gene_id_intersect) == 0 and len(gene_id_combine) == 0:
                    gene_id_combine, gene_id_intersect = test_opposite_strand(alnmt, genes, 2)
 

                # first check the intersection set
                if len(gene_id_intersect) == 1:
                    gene_id = list(gene_id_intersect)[0] + exonintron
                    counts[gene_id] += 1

                elif len(gene_id_intersect) > 1:
                    gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end, RT_barcode) + exonintron
                    counts[gene_id] += 1

                else:
                    # if there no intersection match, then find the union sets
                    if len(gene_id_combine) == 1:
                        gene_id = list(gene_id_combine)[0] + exonintron
                        counts[gene_id] += 1

                    elif len(gene_id_combine) > 1:
                        gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end, RT_barcode) + exonintron
                        counts[gene_id] += 1

                    else:
                        counts["_no_feature"] += 1

    return counts


## Function to assign ambigious reads based on the nearest 3 prime end when it is a shortDT read
def find_nearest_gene(al_end, gene_id_intersect, gene_end, RT_barcode):
    gene_id_end = {}
    for gene in gene_id_intersect:
        if gene in gene_end:
            gene_id_end[gene] = (abs(np.array(list(gene_end[gene])) - al_end)).min()
        else:
            print("****************Found one gene without transcript annotation*****************", "Gene name: ", gene)
    # filter the gene with the least distance. If there are two genes with the least distance, then  "_ambiguous" 
    # would be returned
    gene_end_min = np.min(list(gene_id_end.values()))
    count = 0
    #print ("gene end distance: ", gene_id_end.values())
    for gene in gene_id_end:
        if (gene_id_end[gene] < gene_end_min + 100):
            count += 1
            gene_id = gene
    if count > 1:
        gene_id = "_ambiguous"
    
    if RT_barcode == "randomN":
        gene_id = "_ambiguous"
    
    return(gene_id)


## Function to test opposite strand mapping 
def test_opposite_strand(almnt_tmp, genes_tmp, exongene):
    gene_id_intersect_tmp = set()
    gene_id_combine_tmp = set()
    inter_count_tmp = 0
    
    for cigop in almnt_tmp.cigar:
        if cigop.type != "M":
            continue
        if exongene == 1:
            if cigop.ref_iv.strand == "+":
                cigop.ref_iv.strand = "-"
            else:
                cigop.ref_iv.strand = "+"                
        for iv,val in genes_tmp[cigop.ref_iv].steps():
            gene_id_combine_tmp |= val
            if inter_count_tmp == 0:
                gene_id_intersect_tmp |= val
                inter_count_tmp += 1
            else:
                gene_id_intersect_tmp &= val
                
    return(gene_id_combine_tmp, gene_id_intersect_tmp)


## Function to test 3 prime end flanking regions mapping
def test_flanking_region(almnt_tmp, genes_tmp):
    gene_id_intersect_tmp = set()
    gene_id_combine_tmp = set()
    inter_count_tmp = 0
    
    if almnt_tmp.iv.strand == "+":
        
        if (almnt_tmp.iv.start - 1000) >= 0:
            difference = 1000
        else:
            difference = almnt_tmp.iv.start
     
        for cigop in almnt_tmp.cigar:
            if cigop.type != "M":
                continue

            cigop.ref_iv.start =  cigop.ref_iv.start - difference
            cigop.ref_iv.end =  cigop.ref_iv.end - difference
    
            for iv,val in genes_tmp[cigop.ref_iv].steps():
                gene_id_combine_tmp |= val
                if inter_count_tmp == 0:
                    gene_id_intersect_tmp |= val
                    inter_count_tmp += 1
                else:
                    gene_id_intersect_tmp &= val
            cigop.ref_iv.start = cigop.ref_iv.start + difference
            cigop.ref_iv.end = cigop.ref_iv.end + difference


    else:

        for cigop in almnt_tmp.cigar:
            if cigop.type != "M":
                continue
            cigop.ref_iv.start = cigop.ref_iv.start + 1000
            cigop.ref_iv.end = cigop.ref_iv.end + 1000
            for iv,val in genes_tmp[cigop.ref_iv].steps():
                gene_id_combine_tmp |= val
                if inter_count_tmp == 0:
                    gene_id_intersect_tmp |= val
                    inter_count_tmp += 1
                else:
                    gene_id_intersect_tmp &= val
            cigop.ref_iv.start = cigop.ref_iv.start - 1000
            cigop.ref_iv.end = cigop.ref_iv.end - 1000
    
    return(gene_id_combine_tmp, gene_id_intersect_tmp)

               
if __name__ == "__main__":   
    gtf_file = sys.argv[1]
    input_folder = sys.argv[2]
    output_folder = sys.argv[3]
    sample_ID_file = sys.argv[4]
    core_number = sys.argv[5]
    randomN_barcode_file = sys.argv[6]
    EasySci_count_parallel(gtf_file, input_folder, output_folder, sample_ID_file, core_number, randomN_barcode_file)

