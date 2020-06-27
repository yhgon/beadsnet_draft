%%file get_ids_simple_v12.py
import sys, time, argparse, json, os, shutil, glob

AA_basic ='*ARNDCQEGHILKMFPSTWYV' # 21
AA_ext1  ='*ARNDCQEGHILKMFPSTWYVUO'
AA_ext2  ='*ARNDCQEGHILKMFPSTWYVUOXBZJ_ .'

AAb_to_int = dict((c, i) for i, c in enumerate(AA_basic))
AA1_to_int = dict((c, i) for i, c in enumerate(AA_ext1))
AA2_to_int = dict((c, i) for i, c in enumerate(AA_ext2))

int_to_AAb = dict((i, c) for i, c in enumerate(AA_basic))
int_to_AA1 = dict((i, c) for i, c in enumerate(AA_ext1))
int_to_AA2 = dict((i, c) for i, c in enumerate(AA_ext2))


def seq_to_intlist(seq):
   result =  [AA2_to_int[char] for char in seq]
   return result


def intlist_to_seq(arr):
    result =  [int_to_AA2[code] for code in arr]
    result = ''.join(result)
    return result



def get_fasta_filesize(filename =None):
    file_stats = os.stat(filename)
    file_size_sum = file_stats.st_size
    str_file_size=''
    if file_size_sum >= 1000*1000*1000:
        str_file_size = "{:4.2f}GB".format(  file_size_sum / (1024*1024*1024)   )
    elif file_size_sum < 1000*1000*1000 and file_size_sum >= 1000*1000:
        str_file_size = "{:4.2f}MB".format(  file_size_sum / (1024*1024)   )
    elif file_size_sum < 1000*1000 and file_size_sum > 1000:
        str_file_size = "{:4.2f}KB".format(  file_size_sum / (1024)   )
    elif file_size_sum < 1000 and file_size_sum > 1000:
        str_file_size = "{}Byte".format( file_size_sum   )                 
    output =  str_file_size
    return output

def get_numlines_v12(filename='sample.fasta', tar_seq=0 ):
        
    byte_line = 14+4+6+ 2+1  #27 =id(12),len(3:999),hoffset(14),hlen(4),slen(6),delimiter(4)+newline(1)
    byte_offset= tar_seq*byte_line
    byte_read = byte_line # eliminate last 
    f = open(filename, 'rb')
    f.seek(byte_offset,1)
    line_value = f.read(byte_read) 
    f.close()
    str_value = str(line_value.decode('utf-8').rstrip("\n"))
    line_info = str_value.replace('\n','').replace(" ","").split(",")
    idx       =  tar_seq  # only num_mode=0 
    off_h     =  int(line_info[0])   #2 offset for head
    len_h     =  int(line_info[1]) #3
    off_s     =  off_h + len_h           #
    len_s     =  int(line_info[2])  #4
    return idx, off_h, len_h, off_s, len_s


def header_parse_uniport(header):
    hsplit = header.split('|')
    h_did  = hsplit[0].replace('>','')
    h_uid  = hsplit[1]     
    _des1  = hsplit[2].split('OS=')

    _des2  = _des1[1].split('OX=')
    _des3  = _des2[1].split('GN=')
    _des4  = _des3[1].split('PE=')
    _des5  = _des4[1].split('SV=')

    val_OS = _des2[0]             
    val_OX = _des3[0]             
    val_GN = _des4[0]         
    val_PE = _des5[0]        
    val_SV = _des5[1]

    _desc6   =  _des1[0].split(',')
    _desc7   =  _desc6[0].split(' ')
    val_for  =  _desc7[0]
    val_bac1 =  _desc7[1]
    val_bac2=None
    if len(_desc6)==2:
        val_bac2 = _desc6[1]

    str_val_back = str( _desc6[0] ).strip().replace(str(val_for+' '), '')

    val_bac1
    str1="|{}| {}".format( h_did, h_uid, )
    str2="|{}| {}| {}".format(val_for, str_val_back, val_bac2)
    str3="|{}| {}| {}| {}| {}".format(val_OS, val_OX, val_GN, val_PE, val_SV )  

    #print(str1+str2+str3)
    #print(str_val_back)
    return h_did, h_uid, val_for, val_OS, val_OX, val_GN, val_PE, val_SV, str_val_back, val_bac2
        

def get_header_one_v12(filename='sample.fasta' ,   idx=0,  off_h=0,  len_h=0,  off_s=0,     len_s=0):
    byte_offset = off_h 
    byte_read   = len_h # eliminate last 
    f = open(filename, 'rb')
    f.seek(byte_offset,1)
    header_value = f.read(byte_read) 
    f.close()  
    str_header_value1 = str(header_value.decode('utf-8').rstrip("\n"))
    str_header_value2 = str_header_value1.replace('\n','') 
    return idx, str_header_value2


def get_seq_one_v12(filename='sample.fasta', idx=0,  off_h=0,  len_h=0,  off_s=0, len_s=0):   
    byte_offset = off_s  # off_h + len_h
    byte_read   = len_s  # off_h + len_h + len_s 
    f = open(filename, 'rb')
    f.seek(byte_offset,1)
    byte_value = f.read(byte_read) 
    f.close()
    str_value1 = str(byte_value.decode('utf-8').rstrip("\n"))
    str_value2 = str_value1.replace('\n','') 
    return idx, str_value2


 
def main():
    parser = argparse.ArgumentParser(description="FAI_lite v5 ")
    # Set the default for the dataset argument
    parser.add_argument( "--fasta_dir",   type=str, default='/Protein/UniRef/')
    parser.add_argument( "--idx_dir",   type=str, default='/Protein/UniRef/pkl_uniref50_100k')
    parser.add_argument( "--filename",  type=str, default='uniref50.fasta' )
    parser.add_argument( "--tar_seq",   type=int, default=0 )  

    args = parser.parse_args()
 
    print("============= fai_lite v5. by Hyungon Ryu | NVAITC =============")
    print(" generate line info for target sequence index for DL training ")    
    print("   input  : idx file, target_seq  ")       
    print("   output : line number and length of lines for Seq.\n")          
    print("     args : fasta{} idx{} file{} seq{} ".format(args.fasta_dir, args.idx_dir, args.filename, args.tar_seq  ))
 
    # Create a dictionary of the shell arguments
    import os
 
    if not os.path.exists(args.idx_dir):
        os.makedirs(args.idx_dir)

################################ file name 
    filebody,ext  = os.path.splitext(args.filename)
    fasta_filename = os.path.join(args.fasta_dir, args.filename)    
    idx_filename = os.path.join(args.idx_dir, args.filename +'.idx'+'_v12')    

    
    
    for i in range(130623):                
        idx, off_h, len_h, off_s, len_s = get_numlines_v12  (idx_filename,   tar_seq=i) 
        idx, seq                        = get_seq_one_v12   (fasta_filename, idx, off_h, len_h, off_s, len_s  )   
        idx, header                     = get_header_one_v12   (fasta_filename, idx, off_h, len_h, off_s, len_s  )          
        ids                             = seq_to_intlist(seq )
        
        if i  > 100 :
            return 

        print(header_parse_uniport(header) )

if __name__ == "__main__":
    main()

