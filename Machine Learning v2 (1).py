#!/usr/bin/env python
# coding: utf-8

# # Libraries

# In[2]:


#Import general modules 

from datetime import datetime
import numpy as np
import re
import matplotlib.pyplot as plt 
import pandas as pd


# In[3]:


# import BioPython modules

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import GenBank 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Entrez


# In[4]:


# Import PHASTER modules

import requests
import wget
import webbrowser
import json


# In[5]:


# Import BLAST modules

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline


# In[6]:


# Import  machine learning modules

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
import tensorflow_addons as tfa
import tensorflow_io as tfio

from keras.callbacks import EarlyStopping
from keras.callbacks import ReduceLROnPlateau


# # Import Genome Files

# In[7]:


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
    


# In[8]:


D133 = import_fasta_files('D133.fasta')

E185B = import_fasta_files('E185B.fasta')

E011 = import_fasta_files('E011.fasta')

S4_1 = import_fasta_files('S4-1_genome.fasta')

S8 = import_fasta_files('S8_genome.fasta')


# # Machine Learning

# In[7]:


# printing the version of tensorflow
print("TensorFlow version: ", tf.__version__)


# ## One-Hot encoding

# In[9]:


from sklearn.preprocessing import OneHotEncoder

def one_hot_encoder(sequence):
    encoder = OneHotEncoder()
    x = np.array(list(sequence)).reshape(-1, 1)
    x_onehotencoder = encoder.fit_transform(x).toarray()

    print(x_onehotencoder.shape)

    return x_onehotencoder


# In[10]:


D133_onehot = one_hot_encoder(D133['Code'][0])


# In[11]:


E185B_onehot = one_hot_encoder(E185B['Code'][0])


# In[12]:


E011_onehot = one_hot_encoder(E011['Code'][0])


# In[13]:


S4_1_onehot = one_hot_encoder(S4_1['Code'][0])


# In[14]:


S8_onehot = one_hot_encoder(S8['Code'][0])


# In[15]:


D133_onehot_time = get_ipython().run_line_magic('timeit', '-o D133_onehot')

E185B_onehot_time = get_ipython().run_line_magic('timeit', '-o E185B_onehot')

E011_onehot_time = get_ipython().run_line_magic('timeit', '-o E011_onehot ')

S4_1_onehot_time = get_ipython().run_line_magic('timeit', '-o S4_1_onehot ')

S8_onehot_time = get_ipython().run_line_magic('timeit', '-o S8_onehot ')


# ## Sequential Encoding

# - Encode categorical features between 0 and n_classes-1

# In[152]:


from sklearn import preprocessing
encoder = preprocessing.LabelEncoder()
x = D133['Code'][0]
x1 = np.array(list(x)).reshape(-1, 1)
z = encoder.fit(x1)
y = encoder.transform(x1)

y = (y + 1) / 4
y


# In[153]:


from sklearn import preprocessing

def labelencoder(sequence):
    encoder = preprocessing.LabelEncoder()
    x = sequence
    x1 = np.array(list(x)).reshape(-1, 1) # have to separate out otherwise i get an error for this way
    x2 = encoder.fit(x1)
    x_labelencoder = encoder.transform(x1)
    x_labelencoder = (x_labelencoder + 1) / 4
    
    print(x_labelencoder.shape)

    return x_labelencoder


# In[154]:


D133_labelencoder = labelencoder(D133['Code'][0])


# In[155]:


E185B_labelencoder = labelencoder(E185B['Code'][0])


# In[156]:


E011_labelencoder = labelencoder(E011['Code'][0])


# In[157]:


S4_1_labelencoder = labelencoder(S4_1['Code'][0])


# In[158]:


S8_labelencoder = labelencoder(S8['Code'][0])


# In[159]:


D133_label_time = get_ipython().run_line_magic('timeit', '-o D133_labelencoder')

E185B_label_time = get_ipython().run_line_magic('timeit', '-o E185B_labelencoder')

E011_label_time = get_ipython().run_line_magic('timeit', '-o E011_labelencoder')

S4_1_label_time = get_ipython().run_line_magic('timeit', '-o S4_1_labelencoder')

S8_label_time = get_ipython().run_line_magic('timeit', '-o S8_labelencoder ')


# # Ordinal 

# - Encode categorical features using an ordinal encoding scheme.

# In[103]:


encoder = OrdinalEncoder()
x = D133['Code'][0]
x1 = np.array(list(x)).reshape(-1, 1)
z = encoder.fit(x1)
encoder.categories_
y = encoder.transform(x1)
y


# In[112]:


from sklearn.preprocessing import OrdinalEncoder

def ordinalencoder(sequence):
    encoder = OrdinalEncoder()
    x = sequence
    x1 = np.array(list(x)).reshape(-1, 1) # have to separate out otherwise i get an error for this way
    x2 = encoder.fit(x1)
    x_ordinalencoder = encoder.transform(x1)
    
    print(x_ordinalencoder.shape)

    return x_ordinalencoder


# In[113]:


D133_ordinalencoder = ordinalencoder(D133['Code'][0])


# In[114]:


E185B_ordinalencoder = ordinalencoder(E185B['Code'][0])


# In[115]:


E011_ordinalencoder = ordinalencoder(E011['Code'][0])


# In[116]:


S4_1_ordinalencoder = ordinalencoder(S4_1['Code'][0])


# In[117]:


S8_ordinalencoder = ordinalencoder(S8['Code'][0])


# In[118]:


D133_ordinal_time = get_ipython().run_line_magic('timeit', '-o D133_ordinalencoder')

E185B_ordinal_time = get_ipython().run_line_magic('timeit', '-o E185B_ordinalencoder')

E011_ordinal_time = get_ipython().run_line_magic('timeit', '-o E011_ordinalencoder')

S4_1_ordinal_time = get_ipython().run_line_magic('timeit', '-o S4_1_ordinalencoder')

S8_ordinal_time = get_ipython().run_line_magic('timeit', '-o S8_ordinalencoder ')


# ### K-mer Encoding

# In Bioinformatics, k-mers are subsequences of length contained within a biological sequence. Usually, the term k-mer refers to all of a sequence’s subsequences of length k, such that the sequence AGAT would have four monomers (A, G, A, and T), three 2-mers (AG, GA, AT), two 3-mers (AGA and GAT) and one 4-mer (AGAT) — from k-mer definition. In general, decomposing a sequence into its k-mers fixed-size chunks allows fast and easy string manipulation. This is wisely applied in Natural Language Processing (NLP) bag of words method for ML algorithms. This method will be covered in the next topic.

# In[168]:


def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers


# In[ ]:


phiC2_length = 56538


# In[ ]:


tn6215_length = 


# In[6]:


# Importing the Tn6215 genome into python

# Have to set my email otherwise Entrez does not work
Entrez.email = 'sjjbyrne@hotmail.com'

handle = Entrez.efetch(db="nuccore", id="DQ466086", rettype="gb", retmode="text")
genome = SeqIO.parse(handle, "genbank")
for record in genome:
    print(f'     The record id = \n {record.id} \n \n     The record descriptions = \n {record.description} \n \n     The record = \n {record} \n \n     The record sequence = \n {record.seq}')

phiC2 = record.seq


# In[169]:


D133_kmerencoder = build_kmers(D133['Code'][0], 56538)
D133_kmerencoder


# In[ ]:





# In[ ]:





# ## Visualisation

# In[16]:


# plot table of it
index = ['D133', 'E185B', 'E011', 'S4 1', 'S8']

data = {
    'OneHot Timeit' : [32.4, 30.4, 30.6, 31.2, 34.2],
    'OneHot Error' :[1.7, 0.641, 0.61, 0.713, 3],
    'Label Timeit' : [36, 33, 36, 35.6, 43.3],
    'Label Error' :[2.43, 1.6, 2.04,  3.81, 3.13],
       }

encoding_df = pd.DataFrame(data=data,
                        index=index)
encoding_df


# In[22]:


encoding_df['OneHot Timeit'].sum() / 5 


# In[24]:


encoding_df['Label Timeit'].sum() / 5 


# In[ ]:


encoding_df['Kmer Timeit'].sum() / 5 


# In[29]:


index = ['Mean Speed (ms)']

data = {
    'One Hot': encoding_df['OneHot Timeit'].sum() / 5,
    'Label': encoding_df['Label Timeit'].sum() / 5 
}

df = pd.DataFrame(data=data, 
                 index=index)

df


# In[ ]:





# In[ ]:


# plot table of it
index = ['D133', 'E185B', 'E011', 'S4 1', 'S8']

data = {
    'OneHot Timeit' : [32.4, 30.4, 30.6, 31.2, 34.2],
    'OneHot Error' :[1.7, 0.641, 0.61, 0.713, 3],
    'Label Timeit' : [36, 33, 36, 35.6, 43.3],
    'Label Error' :[2.43, 1.6, 2.04,  3.81, 3.13],
       }

encoding_df = pd.DataFrame(data=data,
                        index=index)
encoding_df


# In[31]:


fig1, axs1 = plt.subplots(1,1,
                         figsize=(10,10))

axs1.errorbar(encoding_df.index.values, 
              encoding_df['OneHot Timeit'], 
              yerr=encoding_df['OneHot Error'], 
              fmt="*",
              color='navy',
              label='One-Hot Encoder')

axs1.errorbar(encoding_df.index.values, 
              encoding_df['Label Timeit'], 
              yerr=encoding_df['Label Error'], 
              fmt=".",
              color='green',
              label='Label Encoder')

axs1.set_ylabel('Time Taken per loop (ns)')
axs1.set_xlabel('Genome Isolate')


plt.legend()

plt.show()


# In[ ]:


Tn6215	PhiC2	
1	5	D133
3	6	E185B
1	5	E011
2	5	S4 1
1	2	S8


# In[ ]:


# plot table of it
index = ['D133', 'E185B', 'E011', 'S4 1', 'S8']

data = {
    'OneHot Timeit' : [32.4, 30.4, 30.6, 31.2, 34.2],
    'OneHot Error' :[1.7, 0.641, 0.61, 0.713, 3],
    'Label Timeit' : [36, 33, 36, 35.6, 43.3],
    'Label Error' :[2.43, 1.6, 2.04,  3.81, 3.13],
       }

encoding_df = pd.DataFrame(data=data,
                        index=index)
encoding_df


# # Split Train etc.

# In[32]:


from sklearn.model_selection import train_test_split


# In[56]:


def test_train_split(x):
    
    x_train, x_test = train_test_split(x, 
                                       test_size=0.05)
    
    print('train =', x_train.shape)
    print('test =', x_test.shape)
    
    return x_train, x_test


# In[45]:


test_train_split(D133_onehot)


# In[46]:


test_train_split(E011_onehot)


# In[47]:


test_train_split(S4_1_onehot)


# In[48]:


test_train_split(E185B_onehot)


# In[49]:


test_train_split(S8_onehot)


# # CNN

# In[51]:


from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, Dense, MaxPooling1D, Flatten
from tensorflow.keras.optimizers import Adam


# In[52]:


model = Sequential()
model.add(Conv1D(64, 3, activation='relu', input_shape=(72, 1)))
model.add(MaxPooling1D(2))
model.add(Conv1D(64, 3, activation='relu'))
model.add(MaxPooling1D(2))
model.add(Flatten())
model.add(Dense(64, activation='relu'))
model.add(Dense(28, activation='softmax'))

model.compile(Adam(lr=.0001), loss='categorical_crossentropy', metrics=['accuracy'])

model.summary()


# In[57]:


S8_train, S8_test = test_train_split(S8_onehot)


# In[61]:


S8_train.shape


# In[58]:


model.fit(S8_train, 
          batch_size=32, 
          epochs=3, 
          validation_split=0.1)


# In[ ]:





# In[ ]:





# In[ ]:




