import os, sys

from itertools import islice

import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torch.utils.data import IterableDataset

import numpy as np 

import argparse, json, yaml

import dataloader * 


parser = argparse.ArgumentParser(description='Speech Synthesis PyTorch Tacotron 2 Training')
args, _ = parser.parse_known_args()   
args.batch_size   = 3
args.fasta_file   = 'complete.nonredundant_protein.451.protein.faa'
args.fasta_dir    = '/content' 
args.idx_dir      = '/content/'
args.idx_file     = 'complete.nonredundant_protein.451.protein.faa.idx_v12'
args.chars_line   = 27
args.on_mem       = False
args.vocab_option = 'ext2'
args.n_seq_total  = 130623 
args.max_len      = 64

print ( json.dumps( vars(args) , sort_keys=False, indent=4 ) ) 


my_data = ProteinDataset(args)
my_collate_fn = Protein_Collate(args.max_len)
data_loader = DataLoader(my_data, num_workers=0, 
                         shuffle=False, 
                         batch_size=args.batch_size, 
                         pin_memory=False,
                         drop_last=True)

import numpy as np
epochs = 2

for i in range(epochs):
    for j, batch in enumerate(data_loader) :
        if j >= 10 :
            break
        print(i, j ,    batch  )
