'''
This script reads in the aligned and filtered sam files, then removes the duplicates based on the cell identity, UMI sequence and mapping locations
'''

import sys
from functools import partial
from multiprocessing import Pool

## This function removes the duplicates based on the cell identity, UMI sequence and mapping locations
def rm_dup_samfile(input_folder, output_folder, sampleID):
    f1 = open(input_folder + sampleID + '.sam')
    f2 = open(output_folder + sampleID + '.sam', 'w')
    pre_location = set()
    pre_barcode_UMI = []    
    pre_dup_num = 0
    cur_dup_num = 0

    line1 = f1.readline()
    
    while(line1):
        
        if (line1[0] == '@'):
            f2.write(line1)
            line1 = f1.readline()
        else:
            name1 = (((line1.split('\t'))[0]).split(','))
            barcode_UMI = name1[0] + name1[1]
            chrom_num = (line1.split('\t'))[2]
            start_site = (line1.split('\t'))[3]            
            current_location = chrom_num + "-" + start_site

            if barcode_UMI == pre_barcode_UMI:    
                dup = False
                if current_location in pre_location:
                    dup = True
                if dup == False:
                    pre_dup_num = cur_dup_num
                    cur_dup_num = 1
                    f2.write(line1)
                    pre_location.add(current_location)
                    line1 = f1.readline()
                else:
                    cur_dup_num += 1
                    line1 = f1.readline()              
            else:
                pre_dup_num = cur_dup_num
                cur_dup_num = 1
                f2.write(line1)
                pre_barcode_UMI = barcode_UMI
                pre_location = set()
                pre_location.add(current_location)
                line1 = f1.readline()
    
    f1.close()
    f2.close()


## This function paralellizes the duplicate removal across PCR batches
def run_paralell(input_folder, sampleID, output_folder, cores):
    
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    p = Pool(processes = int(cores))
    func = partial(rm_dup_samfile, input_folder, output_folder)
    result = p.map(func, sample_list)
    p.close()
    p.join()


if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    cores = sys.argv[4]
    run_paralell(input_folder, sampleID, output_folder, cores)

