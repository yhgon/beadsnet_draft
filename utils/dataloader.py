import os, sys

from itertools import islice

import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torch.utils.data import IterableDataset

import numpy as np 


def load_offsetdata(filename, split=","):
    with open(filename, encoding='utf-8') as f:
        offsetinfo = [line.strip().split(split) for line in f]
    return offsetinfo

class ProteinDataset(Dataset):
    def __init__(self, args ): 

        ########## for index offset ########################
        self.fasta_file = os.path.join(args.fasta_dir, args.fasta_file )
        self.idx_file   = os.path.join(args.idx_dir, args.idx_file )
        self.chars_line = args.chars_line # 27

        self.on_mem     = args.on_mem
        self.offsetinfo = None
        if self.on_mem :
            self.offsetinfo   = load_offsetdata(self.idx_file)     

        ############ for protein vocab ###########################
        self.vocab_option = args.vocab_option  #       
        self.DICT_AA_b    = '*ARNDCQEGHILKMFPSTWYV' # 21
        self.DICT_AA_ext1 = '*ARNDCQEGHILKMFPSTWYVUO' 
        self.DICT_AA_ext2 = '*ARNDCQEGHILKMFPSTWYVUOXBZJ'        
        self.DICT_DSSP    = 'LHBEGITS' # for secondary structure 

        self.AAb_to_int = dict((c, i) for i, c in enumerate(self.DICT_AA_b))  
        self.AA1_to_int = dict((c, i) for i, c in enumerate(self.DICT_AA_ext1))
        self.AA2_to_int = dict((c, i) for i, c in enumerate(self.DICT_AA_ext2))
        self.DSSP_to_int = dict((c, i) for i, c in enumerate(self.DICT_DSSP))        

        self.int_to_AAb = dict((i, c) for i, c in enumerate(self.AA_basic))
        self.int_to_AA1 = dict((i, c) for i, c in enumerate(self.DICT_AA_ext1))
        self.int_to_AA2 = dict((i, c) for i, c in enumerate(self.DICT_AA_ext2))      
        self.int_to_DSSP = dict((i, c) for i, c in enumerate(self.DICT_DSSP))


    def get_offset_disk(self, index ):   
        f = open(self.idx_file, 'rb')
        f.seek(index * self.chars_line,1)
        line_value = f.read(self.chars_line ) 
        f.close()

        str_value = str(line_value.decode('utf-8').rstrip("\n"))
        line_info = str_value.replace('\n','').replace(" ","").split(",")   
        
        off_h     =  int(line_info[0])  
        len_h     =  int(line_info[1])   
        off_s     =  off_h + len_h       
        len_s     =  int(line_info[2])   
    
        return off_h, len_h, off_s, len_s

    def get_offset_mem(self, index ): 
        line_info =  self.offsetinfo[index]
        off_h     =  int(line_info[0])   
        len_h     =  int(line_info[1])  
        off_s     =  off_h + len_h        
        len_s     =  int(line_info[2])    
        return  off_h, len_h, off_s, len_s     

    def get_seq_offset(self, seq_offset, seq_len )        :
        f = open(self.fasta_file, 'rb')
        f.seek(seq_offset,1)
        byte_value = f.read(seq_len) 
        f.close()

        str_value1 = str(byte_value.decode('utf-8').rstrip("\n"))
        str_value2 = str_value1.replace('\n','') 
        return str_value2

    def get_seq_from_index(self, index):
        if self.on_mem :
            off_h, len_h, off_s, len_s = self.get_offset_mem(index)
        else :
            off_h, len_h, off_s, len_s = self.get_offset_disk(index)
        seq = self.get_seq_offset(off_s, len_s)
        ids = self.seq_to_intlist(self.vocab_option, seq)        
        return ids

    def seq_to_intlist(self, vocab_option, seq):
        if vocab_option == 'basic' :
            protein_dict = self.AAb_to_int
        elif vocab_option == 'ext1' :
            protein_dict = self.AA1_to_int
        elif vocab_option == 'ext2' :    
            protein_dict = self.AA2_to_int        

        result =  [protein_dict[char] for char in seq]
        return result

    def intlist_to_seq(self, vocab_option, arr):
        if vocab_option == basic :
            protein_dict = self.int_to_AAb
        elif vocab_option == ext1 :
            protein_dict = self.int_to_AA1
        elif vocab_option == ext2 :    
            protein_dict = self.int_to_AA2        

        result =  [protein_dict[code] for code in arr]
        result = ''.join(result)        
        return result    

    def __len__ (self):
        return int(args.n_seq_total) # it's too large to check length of dataset

    def __getitem__(self, index):
        return self.get_seq_from_index( index ) 

class Protein_Collate():
    def __init__(self, max_len):
        self.max_len = max_len

    def __call__(self, batch):        
        input_lengths, ids_sorted_decreasing = torch.sort(
          torch.LongTensor([len(x[0]) for x in batch]),
          dim=0, descending=True)
        max_input_len = input_lengths[0]

        seq_padded = torch.LongTensor(len(batch), max_input_len)  
        seq_padded.zero_()

        for i in range(len(ids_sorted_decreasing)):
            seq = batch[ids_sorted_decreasing[i]][0]
            seq_padded[i, :seq.size(0)] = seq

        len_x = [x[2] for x in batch]
        len_x = torch.Tensor(len_x)            
        return seq_padded, input_lengths, output_lengths, len_x    


def pad_data(data):  # TODO for collate_fn
    # Find max length of the mini-batch
    max_len = max(list(zip(*data))[0])
    label_list = list(zip(*data))[2]
    txt_list = list(zip(*data))[3]
    padded_tensors = torch.stack([torch.cat((txt, \
            torch.tensor([pad_id] * (max_len - len(txt))).long())) \
            for txt in txt_list])
    return padded_tensors, label_list
                

