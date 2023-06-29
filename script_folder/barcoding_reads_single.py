'''
    this script accept a read1 file (Read 1), a read2 file (I5), a read3 (Read 2) file, a output_file, a ligation barcode list,
    an RT barcode list, then it open the read1, read2 and read3, output file,
    then extract the RT barcode and UMI sequence in the read 1 file and the ligation barcode from the read 2 file, and convert the
    barcode to the real barcode in the barcode list,
    then it attach the barcodes and UMI sequence to the read name of the read3 file (corresponding to Read 2).
'''   

import sys
import gzip
from functools import partial
import pickle
import re
from multiprocessing import Pool
from Bio.Seq import Seq


## This function iterates through the FASTQ files, extracts and corrects the barcodes, trims the barcode, oligo-dT, adapter sequences off and attaches the barcode and UMI sequence to the header of the new FASTQ files
def UMI_attach_read2_barcode_list(sample, input_folder, output_folder, ligation_barcode_list, RT_barcode_list):
      
    # open the read files
    
    # Note: the "R2" fastq file is actually I5, and the "R3" fastq file is actually R2
    R1_input = gzip.open(input_folder + "/" + sample + ".R1.fastq.gz", 'rt', encoding='utf-8')
    R2_input = gzip.open(input_folder + "/" + sample + ".R3.fastq.gz", 'rt', encoding='utf-8')
    I5_input = gzip.open(input_folder + "/" + sample + ".R2.fastq.gz", 'rt', encoding='utf-8')
    
    # open the output files
    R2_output = gzip.open(output_folder + "/" + sample + ".R2.fastq.gz", 'wt',encoding='utf-8')
    
    # Read line 1 of 4 (header) from input files (done before the loop to emulate a do-while loop; we want it to execute at least once. It is done again at the end of the loop)
    R1_in1 = R1_input.readline()
    R2_in1 = R2_input.readline()
    I5_in1 = I5_input.readline()
    
    # Set counts to 0
    total_line = 0
    filtered_line = 0
    
    # Read through files
    while (R1_in1):
        total_line += 1
        
        # Read line 2 of 4 (sequence) from input files
        R1_in2 = R1_input.readline()
        R2_in2 = R2_input.readline() 
        I5_in2 = I5_input.readline()
        
        # Read line 3 of 4 (+ separator) from input files
        R1_in3 = R1_input.readline()
        R2_in3 = R2_input.readline() 
        I5_in3 = I5_input.readline() # not used
        
        # Read line 4 of 4 (quality score) from input files
        R1_in4 = R1_input.readline()
        R2_in4 = R2_input.readline() 
        I5_in4 = I5_input.readline() # not used
        
        # first check if the ligation barcode matches an expected barcode and correct sequence to closest barcode
        target_lig = I5_in2[:10]
        if target_lig in ligation_barcode_list:
            lig_barcode = ligation_barcode_list[target_lig]
            
            # check if the RT barcode matches an expected barcode and correct sequence to closest barcode
            target_RT = R1_in2[8:18]
            if target_RT in RT_barcode_list:

                RT_barcode = RT_barcode_list[target_RT]
                UMI = R1_in2[:8]
               
                # Create header lines for output file including barcodes and UMIs
                R2_out1 = '@' + lig_barcode + RT_barcode + ',' + UMI + ','+ R2_in1[1:]
                    
                # Create and sequence and quality lines for R2              
                R2_out2 = R2_in2
                R2_out4 = R2_in4
                    
                # Trim ends of sequence and quality lines for R2                 
                RT_bar_seq = Seq(target_RT)
                R2_end_seq = str(RT_bar_seq.reverse_complement())
                
                R2_end_loc = re.search(R2_end_seq, R2_out2)
                if R2_end_loc is not None:
                    R2_out2 = R2_out2[:R2_end_loc.start()] + "\n"
                    R2_out4 = R2_out4[:R2_end_loc.start()] + "\n" 
                
                # Pass third separator line to output
                R2_out3 = R2_in3
                
                # Make sure trimmed sequences are still long enough, and write them to output
                if len(R2_out2) > 20:
                    filtered_line += 1
                    
                    R2_output.write(R2_out1)
                    R2_output.write(R2_out2)
                    R2_output.write(R2_out3)
                    R2_output.write(R2_out4)
                
        # Read line 1 of 4 (header) from input files for next cycle
        R1_in1 = R1_input.readline()
        R2_in1 = R2_input.readline()
        I5_in1 = I5_input.readline()
    
    # Close files
    R1_input.close()
    R2_input.close()
    I5_input.close()
    
    R2_output.close()
    
    print("sample name: %s, total line: %f, filtered line: %f, filter rate: %f" 
          %(sample, total_line, filtered_line, float(filtered_line) / float(total_line)))

# this function accept an input folder and a output folder and then generate the output file with the barcode index attached to the FASTQ header
def attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode_file, RT_barcode_file, core):
    
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    ligation barcode file: %s
    RT barcode file: %s
    
    ___________________________________________________________________________
    ''' %(input_folder, sampleID, output_folder, ligation_barcode_file, RT_barcode_file)
    
    print(init_message)
    
    print("Load ligation barcode dictionary...")
    
    # generate the ligation barcode list
    barcodes = open(ligation_barcode_file, "rb")
    ligation_barcode_list = pickle.load(barcodes)
    barcodes.close()
    
    print("Load RT barcode dictionary...")
    
    # generate the RT barcode list:
    barcodes = open(RT_barcode_file, "rb")
    RT_barcode_list = pickle.load(barcodes)
    barcodes.close()
    
    # for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_list to run UMI_attach_read2_barcode_list
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    # parallel for the functions
    p = Pool(processes = int(core))
    func = partial(UMI_attach_read2_barcode_list, input_folder = input_folder, output_folder=output_folder, ligation_barcode_list = ligation_barcode_list, RT_barcode_list=RT_barcode_list)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print(com_message)
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    ligation_barcode_file = sys.argv[4]
    RT_barcode_file = sys.argv[5]
    core=sys.argv[6]
    attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode_file, RT_barcode_file, core)

