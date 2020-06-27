# include line offset for whole seqence
import sys, time, argparse, json, os, shutil, glob

def get_fasta_filesize(filename =None):
    import os

    
    file_stats = os.stat(filename)
    file_size_sum = file_stats.st_size
    if file_size_sum >= 1000*1000*1000:
        str_file_size = "{:4.2f}GB".format(  file_size_sum / (1024*1024*1024)   )
    elif file_size_sum < 1000*1000*1000 and file_size_sum >= 1000*1000:
        str_file_size = "{:4.2f}MB".format(  file_size_sum / (1024*1024)   )
    elif file_size_sum < 1000*1000 and file_size_sum > 1000:
        str_file_size = "{:4.2f}KB".format(  file_size_sum / (1024)   )
    elif file_size_sum < 1000 and file_size_sum > 1000:
        str_file_size = "{}Byte".format( ile_size_sum   )                 
    output =  str_file_size
    return output

def generate_idx_v12(fasta_dir  = '/Protein/UniRef/raw/', 
                                idx_dir  = '/Protein/UniRef/tmp',
                                filename ="uniref50.fasta"):
    """
    Iterates through fasta
    Args:
        fasta (file): fasta_file
    Returns:
        print ( idx, line for header , # of lines for seqs data  ) , eliminate first row
    """

    filebody,ext = os.path.splitext(filename)
    full_filename = os.path.join(fasta_dir, filename)
    idx_filename = os.path.join(idx_dir, filebody+str(ext)+'.idx'+'_v12')
    print("   input  : {} {} ".format( full_filename, get_fasta_filesize(full_filename) )    )
    

    tic=time.time()
    counts =0
    counts_pre = 0
    idx_pre = 0


    offset= 0
    off_h = 0
    len_h = 0
    off_s = 0
    off_h_pre = 0
    len_h_pre = 0
    off_s_pre = 0    
    

    with open(full_filename) as fread:
        with open(idx_filename, "w") as fwrite:
            sequence = ""
            header = None
            for idx,line in enumerate(fread):            
                                
                if line.startswith('>'):     
                    #print(counts, idx, delta , line.replace('\n','')) 
                    off_h = offset      
                    len_h = len(line) 
                    off_s = off_h + len_h 

                    h_stp = off_h_pre
                    h_len = len_h_pre
                    s_stp = off_s_pre
                    s_len = offset - s_stp                    
                    if idx > 0 :
                        #idx_str = "{:12d},{:3d}".format(idx_pre + 1, idx -idx_pre-1 ) # seq start, # of lines for seq 
                        pts = "{:14d},{:4d},{:6d}\n".format(h_stp, h_len,  s_len  ) # s_stp = h_stp + h_len
                        fwrite.write(pts)                                                  

                    idx_pre = idx       #         
                    len_pre = len(line) # save header len  
                    

                                                                
                    off_h_pre = off_h   # previous header's offset start point
                    len_h_pre = len_h   # header's len 
                    off_s_pre = off_s   # previous seq's start point                   
                    counts+=1
                    if ( counts % 1000000 ==0 ):
                        print("log {:8d}00K sequences processed ".format( int(counts/100000 )  ))                    
                offset +=len(line) # for whole loop                      
    toc = time.time()
    dur = toc - tic
    int_dur= int(dur)
    minutes  = int(int_dur / 60 )
    seconds  = int_dur - minutes*60 
    msec = (seconds*1000)
    idx_fsize = get_fasta_filesize(idx_filename)
    new_str = len(str(pts))

    print("   result : {} {} with {:d} sequences dur {:10.4f}sec width {}".format(  idx_filename, idx_fsize, counts, dur, new_str ) )
    return 

def main():
    parser = argparse.ArgumentParser(description="FAI_lite v12 ")
    # Set the default for the dataset argument
    parser.add_argument( "--fasta_dir",   type=str, default='/Protein/UniRef/')
    parser.add_argument( "--idx_dir",   type=str, default='/Protein/UniRef/pkl_uniref50_100k')
    parser.add_argument( "--filename",  type=str, default='uniref50.fasta' )     
    args = parser.parse_args()
 
    print("============= fai_lite v5. by Hyungon Ryu | NVAITC =============")
    print(" generate protein sequence index to preprocess for DL training ")    
    print("   input  : FASTA file with multiseq ")       
    print("   output : idx file for Seq.\n")          
    print("     args : {} {} {} ".format(args.fasta_dir,  args.idx_dir, args.filename ))
 
    # Create a dictionary of the shell arguments
    import os
 
    if not os.path.exists(args.idx_dir):
        os.makedirs(args.idx_dir)
        
    generate_idx_v12(fasta_dir  = args.fasta_dir,  idx_dir  = args.idx_dir ,    filename = args.filename  )
    print("=================================================================")    
        
    

if __name__ == "__main__":
    main()
