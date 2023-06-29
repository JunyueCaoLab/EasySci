'''
This script takes the files from the duplicate removed folder, extracts the single cells and then counts the exon expressions per cell for the paired-end pipeline.
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
    # read in the gtf file, and then construct the genome interval for exons, genes, and gene end dictionary
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

    for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=True ):

        if len(bundle) != 1:
            continue  # Skip multiple alignments
        first_almnt, second_almnt = bundle[0]

        barcode = first_almnt.read.name.split(',')[0]

        if barcode == pre_barcode:
            counts = count_read(first_almnt, second_almnt, genes, counts)
        else:
            if len(counts) == 0:
                counts = count_read(first_almnt, second_almnt, genes, counts)
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
                counts = count_read(first_almnt, second_almnt, genes, counts)
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

    
## Function to assign one read pair to an exon
def count_read(first_almnt, second_almnt, genes, counts):
        
    if (not first_almnt.aligned) or (not second_almnt.aligned):
        counts[ "_unmapped" ] += 2
    elif (first_almnt.iv.chrom not in genes.chrom_vectors) or (second_almnt.iv.chrom not in genes.chrom_vectors):
        counts["_unmapped"] += 2
    else:
        
        if first_almnt.iv.strand == "+":
            first_almnt.iv.strand = "-"
        else:
            first_almnt.iv.strand = "+"        

        gene_id_intersect_1 = set()
        inter_count_1 = 0
        for cigop in first_almnt.cigar:
            if cigop.type != "M":
                continue
            if cigop.ref_iv.strand == "+":
                cigop.ref_iv.strand = "-"
            else:
                cigop.ref_iv.strand = "+"
            for iv,val in genes[cigop.ref_iv].steps():
                if inter_count_1 == 0:
                    gene_id_intersect_1 |= val
                    inter_count_1 += 1
                else:
                    gene_id_intersect_1 &= val

        gene_id_intersect_2 = set()
        inter_count_2 = 0
        for cigop in second_almnt.cigar:
            if cigop.type != "M":
                continue
            for iv,val in genes[cigop.ref_iv].steps():
                if inter_count_2 == 0:
                    gene_id_intersect_2 |= val
                    inter_count_2 += 1
                else:
                    gene_id_intersect_2 &= val


        if len(gene_id_intersect_1) > 0 or len(gene_id_intersect_2) > 0:
            if len(gene_id_intersect_1) == 0:
                counts["_no_feature"] += 1
                genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                if len(np.unique(np.array(genes_2))) == 1:
                    for i in range(0,len(gene_id_intersect_2)):
                        gene_id = list(gene_id_intersect_2)[i]
                        counts[gene_id] += 1/len(gene_id_intersect_2)
                else:
                    counts["_ambiguous"] += 1
            elif len(gene_id_intersect_2) == 0:
                counts["_no_feature"] += 1
                genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                if len(np.unique(np.array(genes_1))) == 1:
                    for i in range(0,len(gene_id_intersect_1)):
                        gene_id = list(gene_id_intersect_1)[i]
                        counts[gene_id] += 1/len(gene_id_intersect_1)
                else:
                    counts["_ambiguous"] += 1
            elif len(gene_id_intersect_1) == 1 and len(gene_id_intersect_2) == 1:
                genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                if genes_1 == genes_2:
                        gene_id = list(gene_id_intersect_1)[0]
                        counts[gene_id] += 1
                        gene_id = list(gene_id_intersect_2)[0]
                        counts[gene_id] += 1
                else:
                    counts["_ambiguous"] += 2
            elif len(gene_id_intersect_1) > 1 and len(gene_id_intersect_2) > 1:
                genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]               
                intersection = list(set(genes_1) & set(genes_2))
                if len(intersection) == 1:
                    gene_id_1 = list(np.array(list(gene_id_intersect_1))[list(np.where(np.array(genes_1) == intersection[0])[0])])
                    gene_id_2 = list(np.array(list(gene_id_intersect_2))[list(np.where(np.array(genes_2) == intersection[0])[0])])
                    for i in range(0,len(gene_id_1)):
                        gene_id = list(gene_id_1)[i]
                        counts[gene_id] += 1/len(gene_id_1)
                    for i in range(0,len(gene_id_2)):
                        gene_id = list(gene_id_2)[i]
                        counts[gene_id] += 1/len(gene_id_2)                
                else:
                    counts["_ambiguous"] += 2
            elif len(gene_id_intersect_1) > 1:
                genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                if(genes_2[0] in genes_1):
                    gene_id = list(gene_id_intersect_2)[0]
                    counts[gene_id] += 1
                    gene_id_1 = list(np.array(list(gene_id_intersect_1))[list(np.where(np.array(genes_1) == genes_2[0])[0])])
                    for i in range(0,len(gene_id_1)):
                        gene_id = list(gene_id_1)[i]
                        counts[gene_id] += 1/len(gene_id_1)
                else:
                    gene_id = list(gene_id_intersect_2)[0]
                    counts[gene_id] += 1
                    counts["_ambiguous"] += 1               
            elif len(gene_id_intersect_2) > 1:
                genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                if(genes_1[0] in genes_2):
                    gene_id = list(gene_id_intersect_1)[0]
                    counts[gene_id] += 1
                    gene_id_2 = list(np.array(list(gene_id_intersect_2))[list(np.where(np.array(genes_2) == genes_1[0])[0])])
                    for i in range(0,len(gene_id_2)):
                        gene_id = list(gene_id_2)[i]
                        counts[gene_id] += 1/len(gene_id_2)                 
                else:
                    gene_id = list(gene_id_intersect_1)[0]
                    counts[gene_id] += 1
                    counts["_ambiguous"] += 1         
        else:
            gene_id_intersect_1, gene_id_intersect_2 = test_opposite_strand(first_almnt, second_almnt, genes, 1)    
            if len(gene_id_intersect_1) > 0 or len(gene_id_intersect_2) > 0:
                if len(gene_id_intersect_1) == 0:
                    counts["_no_feature"] += 1
                    genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                    if len(np.unique(np.array(genes_2))) == 1:
                        for i in range(0,len(gene_id_intersect_2)):
                            gene_id = list(gene_id_intersect_2)[i]
                            counts[gene_id] += 1/len(gene_id_intersect_2)
                    else:
                        counts["_ambiguous"] += 1
                elif len(gene_id_intersect_2) == 0:
                    counts["_no_feature"] += 1
                    genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                    if len(np.unique(np.array(genes_1))) == 1:
                        for i in range(0,len(gene_id_intersect_1)):
                            gene_id = list(gene_id_intersect_1)[i]
                            counts[gene_id] += 1/len(gene_id_intersect_1)
                    else:
                        counts["_ambiguous"] += 1
                elif len(gene_id_intersect_1) == 1 and len(gene_id_intersect_2) == 1:
                    genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                    genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                    if genes_1 == genes_2:
                            gene_id = list(gene_id_intersect_1)[0]
                            counts[gene_id] += 1
                            gene_id = list(gene_id_intersect_2)[0]
                            counts[gene_id] += 1
                    else:
                        counts["_ambiguous"] += 2
                elif len(gene_id_intersect_1) > 1 and len(gene_id_intersect_2) > 1:
                    genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                    genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]               
                    intersection = list(set(genes_1) & set(genes_2))
                    if len(intersection) == 1:
                        gene_id_1 = list(np.array(list(gene_id_intersect_1))[list(np.where(np.array(genes_1) == intersection[0])[0])])
                        gene_id_2 = list(np.array(list(gene_id_intersect_2))[list(np.where(np.array(genes_2) == intersection[0])[0])])
                        for i in range(0,len(gene_id_1)):
                            gene_id = list(gene_id_1)[i]
                            counts[gene_id] += 1/len(gene_id_1)
                        for i in range(0,len(gene_id_2)):
                            gene_id = list(gene_id_2)[i]
                            counts[gene_id] += 1/len(gene_id_2)                
                    else:
                        counts["_ambiguous"] += 2
                elif len(gene_id_intersect_1) > 1:
                    genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                    genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                    if(genes_2[0] in genes_1):
                        gene_id = list(gene_id_intersect_2)[0]
                        counts[gene_id] += 1
                        gene_id_1 = list(np.array(list(gene_id_intersect_1))[list(np.where(np.array(genes_1) == genes_2[0])[0])])
                        for i in range(0,len(gene_id_1)):
                            gene_id = list(gene_id_1)[i]
                            counts[gene_id] += 1/len(gene_id_1)
                    else:
                        gene_id = list(gene_id_intersect_2)[0]
                        counts[gene_id] += 1
                        counts["_ambiguous"] += 1               
                elif len(gene_id_intersect_2) > 1:
                    genes_1 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_1)).split(",")]
                    genes_2 = [l.split('-')[0] for l in ','.join(list(gene_id_intersect_2)).split(",")]
                    if(genes_1[0] in genes_2):
                        gene_id = list(gene_id_intersect_1)[0]
                        counts[gene_id] += 1
                        gene_id_2 = list(np.array(list(gene_id_intersect_2))[list(np.where(np.array(genes_2) == genes_1[0])[0])])
                        for i in range(0,len(gene_id_2)):
                            gene_id = list(gene_id_2)[i]
                            counts[gene_id] += 1/len(gene_id_2)                 
                    else:
                        gene_id = list(gene_id_intersect_1)[0]
                        counts[gene_id] += 1
                        counts["_ambiguous"] += 1 
            else:
                counts["_no_feature"] += 2
    
    return counts


## Function to test opposite strand mapping 
def test_opposite_strand(first_almnt_tmp, second_almnt_tmp, genes_tmp, exongene):
    gene_id_intersect_tmp_1 = set()
    inter_count_tmp_1 = 0
    
    for cigop in first_almnt_tmp.cigar:
        if cigop.type != "M":
            continue
        if exongene == 1:
            if cigop.ref_iv.strand == "+":
                cigop.ref_iv.strand = "-"
            else:
                cigop.ref_iv.strand = "+"                
        for iv,val in genes_tmp[cigop.ref_iv].steps():
            if inter_count_tmp_1 == 0:
                gene_id_intersect_tmp_1 |= val
                inter_count_tmp_1 += 1
            else:
                gene_id_intersect_tmp_1 &= val

    gene_id_intersect_tmp_2 = set()
    inter_count_tmp_2 = 0
    for cigop in second_almnt_tmp.cigar:
        if cigop.type != "M":
            continue
        if exongene == 1:
            if cigop.ref_iv.strand == "+":
                cigop.ref_iv.strand = "-"
            else:
                cigop.ref_iv.strand = "+"
        for iv,val in genes_tmp[cigop.ref_iv].steps():
            if inter_count_tmp_2 == 0:
                gene_id_intersect_tmp_2 |= val
                inter_count_tmp_2 += 1
            else:
                gene_id_intersect_tmp_2 &= val
                
    return(gene_id_intersect_tmp_1, gene_id_intersect_tmp_2)


                  
if __name__ == "__main__":   
    gtf_file = sys.argv[1]
    input_folder = sys.argv[2]
    output_folder = sys.argv[3]
    sample_ID_file = sys.argv[4]
    core_number = sys.argv[5]
    EasySci_count_parallel(gtf_file, input_folder, output_folder, sample_ID_file, core_number)

