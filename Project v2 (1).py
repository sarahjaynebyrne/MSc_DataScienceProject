#!/usr/bin/env python
# coding: utf-8

# # Project

# ### Pip install packages
# 

# In[1]:


# Install the tensorflow addons package,
# which has a nice image rotation function
get_ipython().system('pip install tensorflow-addons')


# In[2]:


# downloading the plugin profile because it wasnt working properly 
get_ipython().system('pip install -U tensorboard_plugin_profile')


# In[3]:


get_ipython().system('pip install keras')


# In[4]:


get_ipython().system('pip install Bio')


# In[5]:


get_ipython().system('pip install wget')


# In[6]:


get_ipython().system('pip install tensorflow_io')


# In[7]:


get_ipython().system('pip install tensorflow')


# ## Libraries

# In[8]:


# Import general modules 

from datetime import datetime
import numpy as np
import re
import matplotlib.pyplot as plt 
import pandas as pd


# In[9]:


# import BioPython modules

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import GenBank 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Entrez


# In[10]:


# Import PHASTER modules

import requests
import wget
import webbrowser
import json


# In[11]:


# Import BLAST modules

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline


# In[12]:


# Import  machine learning modules

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
import tensorflow_addons as tfa
import tensorflow_io as tfio

from keras.callbacks import EarlyStopping
from keras.callbacks import ReduceLROnPlateau


# ## Import the datafiles

# #### Functions

# In[13]:


# define function to import files and make a format that can easily be used

def convert_to_dataframe(file):
    '''
    
    '''
    # 1. pieces is the total length of the list/2
    pieces = len(file)/2
    # 2. split the list into format of [consensus number, code] 
    fasta_file_split = np.array_split(file, pieces)
    # 3. make a dataframe
    df = pd.DataFrame(fasta_file_split, 
                      columns=('Cluster', 'Code'))
    
    return df

def import_fasta_files(file):
    '''
    
    '''
    # create empty list to paste the fasta data into
    fasta_data = []
    with open(file, 'r') as fasta_file:
        fasta_data = fasta_file.read()
        fasta_data = fasta_data.strip().split()
    
    # apply dataframe function to the imported fasta file
    fasta_dataframe = convert_to_dataframe(fasta_data)
    return fasta_dataframe
    


# #### D133

# In[14]:


D133 = import_fasta_files('D133.fasta')
display(D133)


# #### E185B

# In[15]:


E185B = import_fasta_files('E185B.fasta')
display(E185B)


# #### E011

# In[16]:


E011 = import_fasta_files('E011.fasta')
display(E011)


# ### S4 1

# In[17]:


S4_1 = import_fasta_files('S4-1_genome.fasta')
display(S4_1)


# #### S8

# In[18]:


S8 = import_fasta_files('S8_genome.fasta')
display(S8)


# ## PHASTER

# In[ ]:





# ### PHASTER API RESULTS

# RUN API SHIT HERE

# In[19]:


accession_number = [
    'ZZ_b1c43f1a7a',
    'ZZ_415949e066',
    'ZZ_cf441f14ee',
    'ZZ_7b0e466496',
    'ZZ_df930342d1'
]

names = [
    'S8 Genome',
    'E011 Genome',
    'E185B Genome',
    'S4_1 Genome',
    'D133 Genome'
]

phaster_df = pd.DataFrame(list(zip(names, accession_number)),
                         columns=['Genome', 'Accession Number'])

display(phaster_df)

#I renamed the files (details, summary, phaster.fna) with the genome name in front of the file externally


# ## Summary Files

# #### Functions

# In[19]:


def summary_to_dataframe(file):
    '''
    
    '''
    
    # 1. pieces is the total length of the list/2
    pieces = len(file)/16
    
    # 2. split the list into format of [consensus number, code] 
    fasta_file_split = np.array_split(file, pieces)
    
    # 3. make a dataframe
    summary_df = pd.DataFrame(fasta_file_split) 

    # 4. assign columns
    summary_df.columns = summary_df.iloc[0]
    
    #5. drop original first row
    summary_df.drop(0, 
                    axis=0, 
                    inplace=True)
    
    return summary_df

def summaryname_1(summaryname):
    '''
    
    '''
    
    summary = open(summaryname).read().split()
    
    del summary[0:329]
    
    summary.pop(17)
    
    summary_df = summary_to_dataframe(summary)
    
    return summary_df

def summaryname_2(summaryname):
    '''
    
    '''
    summary = open(summaryname).read().split()
    
    del summary[0:336]
    summary.pop(17)
    
    summary_df = summary_to_dataframe(summary)
    return summary_df


# ### Genome Summary Files

# #### D133

# In[20]:


D133_summary = summaryname_2('D133_Genome_summary.txt')
display(D133_summary)


# #### E185B

# In[21]:


E185B_summary = summaryname_2('E185B_Genome_summary.txt')
display(E185B_summary)


# #### E011

# In[22]:


E011_summary = summaryname_2('E011_Genome_summary.txt')
display(E011_summary)


# #### S4 1

# In[23]:


S4_1_summary = summaryname_2('S4_1_Genome_summary.txt')
display(S4_1_summary)


# #### S8

# In[24]:


S8_summary = summaryname_1('S8_Genome_summary.txt')
display(S8_summary)


# # Find Phage Regions

# #### Functions

# In[25]:


def nucleotide_position(df):
    
    x = df.str.replace(':', '-')
    
    x = x.str.split('-', expand=True)
    
    for cols in x.columns:
        if cols == 2:
            x.drop(columns=x.columns[0], axis=1, inplace=True)
            x.rename(columns={1:'Start', 2:'End'}, inplace=True)
       
    x.rename(columns={0:'Start', 1:'End'}, inplace=True)
    
    x[['Start', 'End']] = x[['Start', 'End']].apply(pd.to_numeric)
    
    x = x[['Start', 'End']] - 1
    
    return x

def find_phages(df1, df2): 
    
    result = []
    
    for i, row in df1.iterrows():
        phage_code = df2['Code'][0][df1['Start'][i]:df1['End'][i]]
        number = i
        i += 1
        result.append({'phage_number': number, 
                       'phage_code': phage_code})
    
    df3 = pd.DataFrame(result,
                      columns=['phage_number','phage_code'])
    
    return df3
           


# ### Tn6215

# In[26]:


# Importing the Tn6215 genome into python

# Have to set my email otherwise Entrez does not work
Entrez.email = 'sjjbyrne@hotmail.com'

handle = Entrez.efetch(db="nuccore", id="KC166248", rettype="gb", retmode="text")
genome = SeqIO.parse(handle, "genbank")
for record in genome:
    print(f'     The record id = \n {record.id} \n \n     The record descriptions = \n {record.description} \n \n     The record = \n {record} \n \n     The record sequence = \n {record.seq}')

Tn6215 = record.seq


# ### PhiC2

# In[27]:


# Importing the Tn6215 genome into python

# Have to set my email otherwise Entrez does not work
Entrez.email = 'sjjbyrne@hotmail.com'

handle = Entrez.efetch(db="nuccore", id="DQ466086", rettype="gb", retmode="text")
genome = SeqIO.parse(handle, "genbank")
for record in genome:
    print(f'     The record id = \n {record.id} \n \n     The record descriptions = \n {record.description} \n \n     The record = \n {record} \n \n     The record sequence = \n {record.seq}')

phiC2 = record.seq


# #### D133

# In[28]:


D133_nucleotide_position = nucleotide_position(D133_summary['REGION_POSITION'])
display(D133_nucleotide_position)


# In[29]:


D133_phages = find_phages(D133_nucleotide_position, D133)
display(D133_phages)


# #### E185B

# In[30]:


E185B_nucleotide_position = nucleotide_position(E185B_summary['REGION_POSITION'])
display(E185B_nucleotide_position)


# In[31]:


E185B_phages = find_phages(E185B_nucleotide_position, E185B)
display(E185B_phages)


# #### E011

# In[32]:


E011_nucleotide_position = nucleotide_position(E011_summary['REGION_POSITION'])
display(E011_nucleotide_position)


# In[33]:


E011_phages = find_phages(E011_nucleotide_position, E011)
display(E011_phages)


# #### S4 1

# In[34]:


S4_1_nucleotide_position = nucleotide_position(S4_1_summary['REGION_POSITION'])
display(S4_1_nucleotide_position)


# In[35]:


S4_1_phages = find_phages(S4_1_nucleotide_position, S4_1)
display(S4_1_phages)


# #### S8

# In[36]:


S8_nucleotide_position = nucleotide_position(S8_summary['REGION_POSITION'])
display(S8_nucleotide_position)


# In[37]:


S8_phages = find_phages(S8_nucleotide_position, S8)
display(S8_phages)


# ----

# # BLAST over the Internet

# In[39]:


help(NCBIWWW.qblast)


# In[40]:


def blast_over_the_internet(sequence, entrez_query):
    result_handle = NCBIWWW.qblast(
                    program="blastn",               # the program to search
                    database="nt",                  # Which database to search against 
                    sequence=sequence,              # the sequence to search
                    entrez_query=entrez_query
                )
    
    with open('blast.xml', 'w+') as save_to:
        save_to.write(result_handle.read())
        result_handle.close()
        
    result_handle = open('blast.xml', 'r')
    blast_record = NCBIXML.read(result_handle)

    E_VALUE_THRESH = 0.0001
    ct = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            ct += 1
            if hsp.expect < E_VALUE_THRESH:
                print('\n')
                print('ALIGNMENTS')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('E value:', hsp.expect)
                print(hsp.query[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')
            
    print("\n" + "there are", ct, "sequences in the BLAST output")
    
    return 


# ### E185B

# ##### Tn6215

# In[74]:


#blast1 = blast_over_the_internet(E185B_phages['phage_code'][0])
#blast2 = blast_over_the_internet(E185B_phages['phage_code'][1])
#blast3 = blast_over_the_internet(E185B_phages['phage_code'][2])
#blast4 = blast_over_the_internet(E185B_phages['phage_code'][3])
#blast5 = blast_over_the_internet(E185B_phages['phage_code'][4])
#blast6 = blast_over_the_internet(E185B_phages['phage_code'][5])


# ##### PhiC2

# In[75]:


blast1 = blast_over_the_internet(E185B_phages['phage_code'][0], 'DQ466086')
blast1


# In[76]:


get_ipython().run_line_magic('timeit', 'blast1')


# In[84]:


blast2 = blast_over_the_internet(E185B_phages['phage_code'][1], 'DQ466086')
blast2


# In[85]:


get_ipython().run_line_magic('timeit', 'blast2')


# ----

# # BLAST Locally

# - easier to use the command panel functions to do this because it is super confusing otherwise

# ## Making the database for each genome file

# #### D133

# In[40]:


# make a database for BLAST to search against containing the transposon Tn6215
get_ipython().system('makeblastdb -in D133.fasta -dbtype nucl')


# #### E185B

# In[41]:


# make a database for BLAST to search against containing the transposon Tn6215
get_ipython().system('makeblastdb -in E185B.fasta -dbtype nucl')


# #### E011

# In[42]:


# make a database for BLAST to search against containing the transposon Tn6215
get_ipython().system('makeblastdb -in E011.fasta -dbtype nucl')


# #### S4 1

# In[43]:


# make a database for BLAST to search against containing the transposon Tn6215
get_ipython().system('makeblastdb -in S4-1_genome.fasta -dbtype nucl')


# #### S8

# In[44]:


# make a database for BLAST to search against containing the transposon Tn6215
get_ipython().system('makeblastdb -in S8_genome.fasta -dbtype nucl')


# ## Tn6215 across the whole database

# #### D133

# In[45]:


# querying KC166248 against the WHOLE of the D133 database
get_ipython().system('blastn -task blastn -query KC166248.fasta -db D133.fasta -evalue 1e-20 -num_threads 4 -out D133_blast.txt')


# In[47]:


# viewing the BLAST text
get_ipython().run_line_magic('less', 'D133_blast.txt')


# In[196]:


get_ipython().run_line_magic('timeit', '!blastn -task blastn -query KC166248.fasta -db D133.fasta -evalue 1e-20 -num_threads 4 -out D133_blast.txt')


# #### E185B

# In[48]:


# querying KC166248 against the WHOLE of the E185B database
get_ipython().system('blastn -task blastn -query KC166248.fasta -db E185B.fasta -evalue 1e-20 -num_threads 4 -out E185B_blast.txt')


# In[49]:


# viewing the BLAST text
get_ipython().run_line_magic('less', 'E185B_blast.txt')


# In[50]:


get_ipython().run_line_magic('timeit', '!blastn -task blastn -query KC166248.fasta -db E185B.fasta -evalue 1e-20 -num_threads 4 -out E185B_blast.txt')


# #### E011

# In[51]:


# querying KC166248 against the WHOLE of the D133 database
get_ipython().system('blastn -task blastn -query KC166248.fasta -db E011.fasta -evalue 1e-20 -num_threads 4 -out E011_blast.txt')


# In[52]:


# viewing the BLAST text
get_ipython().run_line_magic('less', 'E185B_blast.txt')


# In[53]:


get_ipython().run_line_magic('timeit', '!blastn -task blastn -query KC166248.fasta -db E011.fasta -evalue 1e-20 -num_threads 4 -out E011_blast.txt')


# #### S4 1

# In[54]:


# querying KC166248 against the WHOLE of the D133 database
get_ipython().system('blastn -task blastn -query KC166248.fasta -db S4-1_genome.fasta -evalue 1e-20 -num_threads 4 -out S4-1_genome_blast.txt')


# In[55]:


# viewing the BLAST text
get_ipython().run_line_magic('less', 'E185B_blast.txt')


# In[56]:


get_ipython().run_line_magic('timeit', '!blastn -task blastn -query KC166248.fasta -db S4-1_genome.fasta -evalue 1e-20 -num_threads 4 -out S4-1_genome_blast.txt')


# #### S8

# In[57]:


# querying KC166248 against the WHOLE of the D133 database
get_ipython().system('blastn -task blastn -query KC166248.fasta -db S8_genome.fasta -evalue 1e-20 -num_threads 4 -out S8_genome_blast.txt')


# In[58]:


# viewing the BLAST text
get_ipython().run_line_magic('less', 'E185B_blast.txt')


# In[59]:


get_ipython().run_line_magic('timeit', '!blastn -task blastn -query KC166248.fasta -db S8_genome.fasta -evalue 1e-20 -num_threads 4 -out S8_genome_blast.txt')


# ----

# ## Visualisation

# In[163]:


from matplotlib import image


# In[166]:


D133_img = image.imread("D133_dnaplotter.png")
E185B_img = image.imread("E185B_dnaplotter.png")
E011_img = image.imread("E011_dnaplotter.png")
S4_1_img = image.imread("S4_1_dnaplotter.png")
S8_img = image.imread("S8_dnaplotter.png")


# In[194]:


fig1,(axs1, axs2, axs3, axs4, axs5) = plt.subplots(1, 5,
                                                    figsize=(30, 5))

fig1.suptitle('Figure 1 - TITLE', fontname="Times New Roman", size=25, fontweight="bold")

# plotting the images
axs1.imshow(D133_img)
axs2.imshow(E185B_img)
axs3.imshow(E011_img)
axs4.imshow(S4_1_img)
axs5.imshow(S8_img)

# set axes titles
axs1.set_title('A', fontname="Times New Roman", size=25, fontweight="bold")
axs2.set_title('B', fontname="Times New Roman", size=25, fontweight="bold")
axs3.set_title('C', fontname="Times New Roman", size=25, fontweight="bold")
axs4.set_title('D', fontname="Times New Roman", size=25, fontweight="bold")
axs5.set_title('E', fontname="Times New Roman", size=25, fontweight="bold")

# remove x-axes from subplots
axs1.get_xaxis().set_visible(False)
axs2.get_xaxis().set_visible(False)
axs3.get_xaxis().set_visible(False)
axs4.get_xaxis().set_visible(False)
axs5.get_xaxis().set_visible(False)

# remove y-axes from subplots
axs1.get_yaxis().set_visible(False)
axs2.get_yaxis().set_visible(False)
axs3.get_yaxis().set_visible(False)
axs4.get_yaxis().set_visible(False)
axs5.get_yaxis().set_visible(False)

# set tight layout
fig1.tight_layout()

# show the images
plt.show()


# In[195]:


# save image

fig1.savefig('figure1.png')


# ---

# # Machine Learning

# In[24]:


# printing the version of tensorflow
print("TensorFlow version: ", tf.__version__)


# ### 1. One-Hot Encoding

# In[60]:


from sklearn.preprocessing import OneHotEncoder

def one_hot_encoder(sequence):
    encoder = OneHotEncoder()
    x = np.array(list(sequence)).reshape(-1, 1)
    x_onehotencoder = encoder.fit_transform(x).toarray()

    print(x_onehotencoder.shape)

    return x_onehotencoder


# ### Tn6215

# In[61]:


# Tn6215
Tn6215_onehot = one_hot_encoder(Tn6215)


# #### D133

# In[64]:


D133_onehot = one_hot_encoder(D133['Code'][0])


# #### E185B

# In[65]:


E185B_onehot = one_hot_encoder(E185B['Code'][0])


# #### E011

# In[66]:


E011_onehot = one_hot_encoder(E011['Code'][0])


# #### S4 1

# In[67]:


S4_1_onehot = one_hot_encoder(S4_1['Code'][0])


# #### S8

# In[68]:


S8_onehot = one_hot_encoder(S8['Code'][0])


# ---

# ### Random Visualisation

# In[123]:


details = {'frequency_of_Tn6215_regions': [1, 3, 1, 2, 1],
 'frequency_of_phiC2_regions': [5, 6, 5, 5, 2],
 'frequency_of_other_phage_regions': [3, 0, 5, 2, 3]}

df45 = pd.DataFrame(details,
                   index=['D133', 'E185B', 'E011', 'S4_1', 'S8'])
df45


# In[132]:


df45.corr()


# In[162]:


fig1, axs1 = plt.subplots(figsize=(10, 15))

labels = df45.index.values
x = np.arange(len(labels)) 

bar_width = 0.2

axs1.bar(
         x-0.3,
         df45['frequency_of_Tn6215_regions'],
         
         label='Tn6215 regions',
         width=bar_width
         )
axs1.bar(
         x,
         df45['frequency_of_phiC2_regions'],
         label='PhiC2 regions',
         width=bar_width
         )
axs1.bar(
         x+0.3,
         df45['frequency_of_other_phage_regions'],
         label='Misc phage regions',
         width=bar_width
         )

axs1.set_ylabel('Frequency')
axs1.set_xlabel('Genomes')

axs1.set_xticks(range(0, 5, 1))
axs1.set_xticklabels(df45.index)

plt.legend()

plt.show()


# ---

# In[28]:


from Bio import AlignIO

from Bio.Align.Applications import ClustalwCommandline


# In[6]:


get_ipython().system('pip install swalign')


# In[7]:


import swalign


# In[29]:


dna_string = x
reference_string = y
match_score = 2
mismatch_score = -1
matrix = swalign.NucleotideScoringMatrix(match_score, mismatch_score)
lalignment_object = swalign.LocalAlignment(matrix)
alignment_object = lalignment_object.align(dna_string, reference_string)
alignment_object.dump()


# ---

# In[30]:


from skbio.alignment import global_pairwise_align_nucleotide


# In[ ]:


x = E185B_phages['phage_code'][5]
y = Tn6215

s1 = x
s2 = y
r = global_pairwise_align_nucleotide(s1, s2)
print(r.to_fasta())


# ---

# In[ ]:


# build a alignment neural network 

# input = two sequences - can be strings or other arrays of data

# out put is to produce an alignment which pairs up elements of the sequence 


# In[108]:


from Bio import pairwise2
from Bio.pairwise2 import format_alignment
sequence_1 = E185B_phages['phage_code'][5]
sequence_2 = Tn6215
print("Needleman-Wunsch Global Algorithm")
alignments = pairwise2.align.globalxx(sequence_1, sequence_2)
for item in alignments:
    print(format_alignment(*item))
print("Smith-Waterman Local Algorithm")
alignments = pairwise2.align.localxx(sequence_1, sequence_2)
for item in alignments:
    print(format_alignment(*item))


# In[ ]:





# ### DNA Sequence Similarity using different ML Algorithms

# In[34]:


get_ipython().system('pip install Squiggle')


# In[63]:


y = E185B_onehot_phage1

x = Tn6215


# In[ ]:


def alignment(x, y):
    
    


# In[65]:





# In[67]:





# In[104]:




#### The Needleman-Wunsch Algorithm

def needleman_wunsch(x, y, match = 1, mismatch = 1, gap = 1):
    nx = len(x)
    ny = len(y)
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:,0] = np.linspace(0, -nx * gap, nx + 1)
    F[0,:] = np.linspace(0, -ny * gap, ny + 1)
    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3
    P[0,:] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if x[i] == y[j]:
                t[0] = F[i,j] + match
            else:
                t[0] = F[i,j] - mismatch
            t[1] = F[i,j+1] - gap
            t[2] = F[i+1,j] - gap
            tmax = np.max(t)
            F[i+1,j+1] = tmax
            if t[0] == tmax:
                P[i+1,j+1] += 2
            if t[1] == tmax:
                P[i+1,j+1] += 3
            if t[2] == tmax:
                P[i+1,j+1] += 4
    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rx.append(x[i-1])
            ry.append(y[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rx.append(x[i-1])
            ry.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j-1])
            j -= 1
    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return '\n'.join([rx, ry])


# In[106]:


x = E185B_phages['phage_code'][5]
y = Tn6215

needleman_wunsch(x, y, match = 1, mismatch = 1, gap = 1)


# In[ ]:


### 2. 

# Sequential Model Definition
def create_model_unoptimised():
    model = Sequential([
    layers.Input(shape=IMG_SIZE+(3,), 
                 name='Input'),
    layers.Conv2D(16, 5, 
                  padding='same', 
                  activation='relu',
                  name='Conv_1'),
    layers.MaxPooling2D(name='Pool_1'),
    layers.Conv2D(32, 4, 
                  padding='same', 
                  activation='relu',
                  name='Conv_2'),
    layers.MaxPooling2D(name='Pool_2'),
    layers.Conv2D(128, 3, 
                  padding='same', 
                  activation='relu',
                  name='Conv_3'),
    layers.MaxPooling2D(name='Pool_3'),
    layers.Flatten(name='Flatten'),
    layers.Dense(512, 
                 activation='relu', 
                 name='dense_1'),
    layers.Dense(1, 
                 activation='sigmoid', 
                 name='Output')
    ], name='CNN')

    return model


# In[ ]:


# Create a version of the model and print the summary
model_unoptimised = create_model_unoptimised()

model_unoptimised.summary()


# In[ ]:


# Adam optimiser
opt = tf.keras.optimizers.Adam(learning_rate=1e-4)
# Binary classification loss
loss_obj = tf.keras.losses.BinaryCrossentropy(from_logits=False)
# Accuracy metric
metrics = ['accuracy']

# Compile model
model_unoptimised.compile(optimizer=opt,
                          loss=loss_obj,
                          metrics=metrics)


### Use the .evaluate method to check the initial loss and accuracy of the model
### on the test dataset
model_unoptimised.evaluate(test_ds)

# Create a TensorBoard callback
logs = "logs/" + datetime.now().strftime("%Y%m%d-%H%M%S")

tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=logs,
                                                      histogram_freq=1,
                                                      profile_batch='29, 291')



# In[ ]:


# Train the Model

# saving the weights as history_unoptimised 
history_unoptimised = model_unoptimised.fit(train_ds,
                                            validation_data=validation_ds,
                                            epochs=EPOCHS,
                                            batch_size=BATCH_SIZE,
                                            callbacks=[tensorboard_callback])


# a dictionary of the losses and metrics at each epoch
history_unoptimised.history 

# check the values available with:
history_unoptimised.history.keys()
history_unoptimised.history.values()

# storing history_unoptimised as a dataframe to download
history_unoptimised_df = pd.DataFrame(data=history_unoptimised.history.values(),
                                      index=history_unoptimised.history.keys())

# downloading the dataframe 
history_unoptimised_df.to_csv('history_unoptimised_df.csv')

# save weights 
model_unoptimised.save_weights('weights.h5')


# In[ ]:


### plot the loss function and metric of model_unoptimised 

# define figure space and axis
fig1, (axs1, axs2) = plt.subplots(1, 2, 
                                  figsize=(10, 6))

# plot loss data on axs1 
axs1.plot(history_unoptimised.history['loss'],
          color='hotpink',
          label='Loss')
axs1.plot(history_unoptimised.history['val_loss'],
          color='thistle',
          label='Validation Loss')

# plot accuracy data on axs2
axs2.plot(history_unoptimised.history['accuracy'],
          color='magenta',
          label='Accuracy')
axs2.plot(history_unoptimised.history['val_accuracy'],
          color='darkmagenta',
          label='Validation Accuracy')

# x-axis 
axs1.set_xlabel('Epoch (#)')
axs2.set_xlabel('Epoch (#)')
axs1.set_xticks(np.arange(0, 25, 5))
axs2.set_xticks(np.arange(0, 25, 5))

# y-axis
axs1.set_ylabel('Loss')
axs2.set_ylabel('Accuracy')

# legend
axs1.legend()
axs2.legend()

# grid
axs1.grid(0.01)
axs2.grid(0.01)

# add overall title to figure1
fig1.suptitle('Figure 1 - The loss and accuracy of Unoptimised Convolutional Model')

# show the figure
plt.show()


# In[ ]:


# plot tensorflow board 

# Load the TensorBoard notebook extension.
get_ipython().run_line_magic('load_ext', 'tensorboard')


# In[ ]:


# Launch TensorBoard and navigate to the Profile tab to view performance profile
get_ipython().run_line_magic('tensorboard', '--logdir=logs')


# In[ ]:





# In[ ]:





# In[ ]:





# In[32]:




models = {}

# Logistic Regression
from sklearn.linear_model import LogisticRegression
models['Logistic Regression'] = LogisticRegression()

# Support Vector Machines
from sklearn.svm import LinearSVC
models['Support Vector Machines'] = LinearSVC()

# Decision Trees
from sklearn.tree import DecisionTreeClassifier
models['Decision Trees'] = DecisionTreeClassifier()

# Random Forest
from sklearn.ensemble import RandomForestClassifier
models['Random Forest'] = RandomForestClassifier()

# Naive Bayes
from sklearn.naive_bayes import GaussianNB
models['Naive Bayes'] = GaussianNB()

# K-Nearest Neighbors
from sklearn.neighbors import KNeighborsClassifier
models['K-Nearest Neighbor'] = KNeighborsClassifier()


# In[33]:


from sklearn.metrics import accuracy_score, precision_score, recall_score

accuracy, precision, recall = {}, {}, {}

for key in models.keys():
    
    # Fit the classifier model
    models[key].fit(x, y)
    
    # Prediction 
    predictions = models[key].predict(x)
    
    # Calculate Accuracy, Precision and Recall Metrics
    accuracy[key] = accuracy_score(predictions, y)
    precision[key] = precision_score(predictions, y)
    recall[key] = recall_score(predictions, y)

