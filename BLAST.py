#!/usr/bin/env python
# coding: utf-8

# # Libraries

# In[20]:


import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import GenBank 
from Bio.SeqIO.FastaIO import SimpleFastaParser

# to use BLAST library
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# In[2]:


import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt


# In[3]:


help(NCBIWWW.qblast) 


# - using nucleotide blast (blastn)
# - tN6215 ACCESSION NUMBER = KC166248.1

# In[4]:


def import_fasta_files(file):
    '''
    
    '''
    # create empty list to paste the fasta data into
    fasta_data = []
    with open(file, 'r') as fasta_file:
        fasta_data = fasta_file.read()
        fasta_data = fasta_data.strip().split()

    return fasta_data


# In[13]:


# D133

D133 = import_fasta_files('D133.fasta')


# In[7]:


result_handle = NCBIWWW.qblast("blastn", "nt", D133[11]) 
result_handle 


# In[8]:


# opening accession number 


# In[9]:


home=D133[11] 


# In[29]:


def blastn_search(sequence):
    result = NCBIWWW.qblast(
                    program="blastn",               # the program to search
                    database="nt",                  # Which database to search against 
                    sequence=sequence,              # the sequence to search
                    alignments=100,                 # number of alignments to show
                    hitlist_size=30,                # Number of hits to return. Default 50 
    )
    
    with open('results.xml', 'w') as save_file: 
        blast_results = NCBIXML.read(result)
        
    return blast_results

#save_file.write(blast_results)


# In[ ]:


blastn_search(home)


# In[25]:


E_VALUE_THRESH = 1e-20 
for record in NCBIXML.parse(open("results.xml")): 
    if record.alignments: 
        print("\n") 
        print("query: %s" % record.query[:100]) 
    for align in record.alignments: 
        for hsp in align.hsps: 
            if hsp.expect < E_VALUE_THRESH: 
                print("match: %s " % align.title[:100])


# In[ ]:





# In[ ]:




