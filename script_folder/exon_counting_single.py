'''
This script takes the files from the duplicate removed folder, extracts the single cells and then counts the exon expressions per cell for the single-end pipeline.
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
def EasySci_count_parallel(gtf_file, input_folder, output_folder, sample_ID_file, core_number):
    # read in the gtf file, and then construct the genome interval for exons
    gtf_file = HTSeq.GFF_Reader(gtf_file, end_included=True)
    gene_annotat_file = output_folder + "/gene_name_annotate.csv"  
    gene_annotat = open(gene_annotat_file, "w")
    
    genes = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    gene_n = 0
    gene_count = 0
    print("Start generating exon genomic arrays....")
  
    for feature in gtf_file:
        if feature.type == "gene":
            gene_n +=1
            genes[ feature.iv ] += feature.attr["gene_id"]
            gene_count += 1

            message = (feature.attr["gene_id"] + "," + feature.attr["gene_type"] + "," 
                        + feature.attr["gene_name"] + "," + str(gene_count) + "\n")

            gene_annotat.write(message)


    print("Detected exon number: ", gene_n)

    gene_annotat.close()

    gene_annotat = pd.read_csv(gene_annotat_file, header=None)
    gene_annotat.index =  gene_annotat[0]

    sample_ID = list(pd.read_csv(sample_ID_file, header=None)[0])
    
    # parallele for the functions
    p = Pool(processes = int(core_number))
    func = partial(EasySci_count, input_folder, output_folder, genes, gene_annotat)
    result = p.map(func, sample_ID)
    p.close()
    p.join()

    print("All analysis done~")
    
    
## Function to separate, count and save reads from different cells
def EasySci_count(input_folder, output_folder, genes, gene_annotat, sample):    
    
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
            counts = count_read(alnmt, genes, counts)
        else:
            if len(counts) == 0:
                counts = count_read(alnmt, genes, counts)
                pre_barcode = barcode
            else:
                unmatched_reads = 0
                ambiguous_reads = 0
                all_reads = 0
                for gene in counts:
                    if (gene in ["_unmapped", "_ambiguous", "_no_feature"]):
                        if gene in ["_no_feature"]:
                            unmatched_reads += counts[gene]
                        else:
                            ambiguous_reads += counts[gene]
                    else:
                        line = str(gene_annotat.loc[gene,3]) + "," + str(cell_ID) + "," + str(counts[gene]) + "\n"
                        count_output.write(line)
                        all_reads += counts[gene]

                unmatched_rate = unmatched_reads / (unmatched_reads + all_reads + ambiguous_reads)
                line_cell = sample + "." + pre_barcode + "," + str(unmatched_rate) + "\n"
                cell_output.write(line_cell)

                counts = collections.Counter()
                cell_ID += 1         
                counts = count_read(alnmt, genes, counts)
                pre_barcode = barcode

    unmatched_reads = 0
    ambiguous_reads = 0
    all_reads = 0
    for gene in counts:
        if (gene in ["_unmapped", "_ambiguous", "_no_feature"]):
            if gene in ["_no_feature"]:
                unmatched_reads += counts[gene]
            else:
                ambiguous_reads += counts[gene]
        else:
            line = str(gene_annotat.loc[gene,3]) + "," + str(cell_ID) + "," + str(counts[gene]) + "\n"
            count_output.write(line)
            all_reads += counts[gene]

    unmatched_rate = unmatched_reads / (unmatched_reads + all_reads + ambiguous_reads)
    line_cell = sample + "." + pre_barcode + "," + str(unmatched_rate) + "\n"
    cell_output.write(line_cell)

    count_output.close()
    cell_output.close()
    
    return 0  

    
## Function to assign one read to an exon, if the read maps to multiple exons within one gene, the read in split between exons, if it maps to multiple genes, the read is counted as ambiguous.
def count_read(alnmt, genes, counts):
        
    if not alnmt.aligned:
        counts[ "_unmapped" ] += 1
    elif alnmt.iv.chrom not in genes.chrom_vectors:
        counts["_unmapped"] += 1
    else:
        
        # First check the intersecting exons (only total overlap with exons are counted)
        gene_id_intersect = set()
        inter_count = 0
        for cigop in alnmt.cigar:
            if cigop.type != "M":
                continue
            for iv,val in genes[cigop.ref_iv].steps():
                if inter_count == 0:
                    gene_id_intersect |= val
                    inter_count += 1
                else:
                    gene_id_intersect &= val


        if len(gene_id_intersect) > 0:
            genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect)).split(",")]
            if len(np.unique(np.array(genes_2))) == 1:
                for i in range(0,len(gene_id_intersect)):
                    gene_id = list(gene_id_intersect)[i]
                    counts[gene_id] += 1/len(gene_id_intersect)
            else:
                counts["_ambiguous"] += 1
        else:
            gene_id_intersect = test_opposite_strand(alnmt, genes)    
            if len(gene_id_intersect) > 0:
                genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect)).split(",")]
                if len(np.unique(np.array(genes_2))) == 1:
                    for i in range(0,len(gene_id_intersect)):
                        gene_id = list(gene_id_intersect)[i]
                        counts[gene_id] += 1/len(gene_id_intersect)
                else:
                    counts["_ambiguous"] += 1
            else:
                counts["_no_feature"] += 1
    
    return counts


## Function to test opposite strand mapping 
def test_opposite_strand(almnt_tmp, genes_tmp):
    gene_id_intersect_tmp = set()
    inter_count_tmp = 0
    
    for cigop in almnt_tmp.cigar:
        if cigop.type != "M":
            continue
        if cigop.ref_iv.strand == "+":
            cigop.ref_iv.strand = "-"
        else:
            cigop.ref_iv.strand = "+"                
        for iv,val in genes_tmp[cigop.ref_iv].steps():
            if inter_count_tmp == 0:
                gene_id_intersect_tmp |= val
                inter_count_tmp += 1
            else:
                gene_id_intersect_tmp &= val
                
    return(gene_id_intersect_tmp)

               
if __name__ == "__main__":   
    gtf_file = sys.argv[1]
    input_folder = sys.argv[2]
    output_folder = sys.argv[3]
    sample_ID_file = sys.argv[4]
    core_number = sys.argv[5]
    EasySci_count_parallel(gtf_file, input_folder, output_folder, sample_ID_file, core_number)

