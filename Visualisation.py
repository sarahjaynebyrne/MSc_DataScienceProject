#!/usr/bin/env python
# coding: utf-8

# # Libraries

# In[1]:


import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt


# In[3]:


import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import GenBank 
from Bio.SeqIO.FastaIO import SimpleFastaParser


# In[4]:


genome_summary = pd.read_excel('Genome summary for SJByrne.xlsx', header=3)
display(genome_summary)


# # Importing the datafiles

# In[5]:


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
    


# In[6]:


# D133
D133 = import_fasta_files('D133.fasta')
display(D133)


# In[7]:


# E011
E011 = import_fasta_files('E011.fasta')
display(E011)


# In[8]:


#E185B
E185B = import_fasta_files('E185B.fasta')
display(E185B)


# In[9]:


# S4_1
S4_1 = import_fasta_files('S4-1_genome.fasta')
display(S4_1)


# In[10]:


# S8
S8 = import_fasta_files('S8_genome.fasta')
display(S8)


# # Performing sequence features on each data file

# ### Functions for Sequence Features

# In[11]:


def sequence_features(df):
    '''
    This function:
    
    Parameters:
    -----------
    
    '''
    # 1. count the number of bases in the sequence consensus
    df['Length of Code'] = df['Code'].str.count("")
    
    # 2. count the total number of bases in each species
    total_sequence_length = sum(df['Length of Code'])
    print(f'The length of this species is {total_sequence_length}')
    
    # 3. check for sequence composisiton of each nucleotide base
    df['A'] = df['Code'].str.count('A')
    df['T'] = df['Code'].str.count('T')
    df['G'] = df['Code'].str.count('G')
    df['C'] = df['Code'].str.count('C')
    
    df['GC_content'] = ((df['G'] + df['C']) /  df['Length of Code']) * 100
    df['AG_content'] = ((df['A'] + df['T']) /  df['Length of Code']) * 100
    
    # 4. show the sequence compoisition total 
    
    # 5. check for total sequence composisiton and total GC content in each species
    df_GC = df['G'] + df['C']
    GC = (sum(df_GC) / total_sequence_length)*100
    print(f'The Total GC % content of this species is {GC}')
           
    # 6. remove the Code column because un-needed
    df = df.drop('Code', axis=1)
    return df


# In[12]:


def sequence_feature_dataframe(df):
    '''
    This function:
    
    Parameters:
    -----------
    '''
    A = df['Code'].str.count('A')
    T = df['Code'].str.count('T')
    G = df['Code'].str.count('G')
    C = df['Code'].str.count('C')
    
    seqfeatures_data = np.array([(sum(A)), (sum(T)), (sum(G)), (sum(C))])
    
    seqfeatures_df = pd.DataFrame(data=seqfeatures_data, 
                              columns=['Total'], 
                              index=['A', 'T', 'G', 'C'])
    return seqfeatures_df


# In[64]:


def length_of_code_dataframe(df):
    '''
    This function:
    
    Parameters:
    -----------
    '''
    df['Length of Code'] = df['Code'].str.count("")
    
    total_sequence_length = sum(df['Length of Code'])
    
    total_sequence_length_df = pd.DataFrame(data=total_sequence_length, 
                                            columns=['Species'], 
                                            index=['Nucleotide Base Length'])
    
    return total_sequence_length_df


# ### Sequence Features on the Genome files

# In[14]:


sequence_features_D133 = sequence_features(D133)
display(sequence_features_D133)


# In[15]:


sequence_features_E011 = sequence_features(E011)
display(sequence_features_E011)


# In[16]:


sequence_features_E185B = sequence_features(E185B)
display(sequence_features_E185B)


# In[17]:


sequence_features_S4_1 = sequence_features(S4_1)
display(sequence_features_S4_1)


# In[18]:


sequence_features_S8 = sequence_features(S8)
display(sequence_features_S8)


# # Visualisation

# ### Functions

# In[19]:


def plot_sequence_features(df):
    '''
    This function:
    
    Parameters:
    -----------
    '''
    # plot the sequence composition as a bar chart
    fig2, (axs2a, axs2b, axs2c, axs2d) = plt.subplots(1, 4, 
                                                  figsize=(20, 8))
    
    fig2.suptitle('DNA Nucleotide composition in each cluster',
             fontsize=15)
    colors = ['aquamarine', 'mediumaquamarine', 'mediumturquoise',  'lightseagreen', 'teal', 'darkslategrey']
    
    # the data
    axs2a.bar(df.index.values, 
             df['A'],
             color=colors)
    axs2b.bar(df.index.values, 
             df['T'],
             color=colors)
    axs2c.bar(df.index.values, 
             df['G'],
             color=colors)
    axs2d.bar(df.index.values, 
             df['C'],
             color=colors)
    
    # individual titles
    axs2a.set_title('A')
    axs2b.set_title('T')
    axs2c.set_title('G')
    axs2d.set_title('C')


    # Set x axis label
    fig2.text(0.5, 0.04, 
              'Nucleotide Base', 
              ha='center', 
              va='center',
              fontsize=12)

    # Set y axis label
    fig2.text(0.06, 0.5, 
              'Amount', 
              ha='center', 
              va='center', 
              rotation='vertical',
              fontsize=12)
    
    # show the figure
    plt.show()


# In[20]:


def length_of_cluster(df):

    fig5, axs5 = plt.subplots(1, 1,
                           figsize=(10,5))

    fig5.suptitle('DNA Length of code in each cluster',
                 fontsize=15)
    
    colors = ['lightseagreen', 'mediumturquoise', 'mediumaquamarine', 'aquamarine']

    axs5.bar(df.index.values,
              df['Length of Code'],
              color=colors)
    
    axs5.set_xlabel('Cluster Number')
    
    axs5.set_ylabel('Nucleotide Length')

    plt.show()


# ### Figures

# In[21]:


# set up the figure 
fig1, (axs1a, axs1b, axs1c, axs1d, axs1e) = plt.subplots(1, 5, 
                                                  figsize=(20, 8))
fig1.suptitle('Figure 1 - DNA Nucleotide composition',
             fontsize=15)

colors = ['lightseagreen', 'mediumturquoise', 'mediumaquamarine', 'aquamarine']

# the data
axs1a.bar(sequence_feature_dataframe(D133).index.values, 
         sequence_feature_dataframe(D133)['Total'],
         color=colors)
axs1b.bar(sequence_feature_dataframe(E011).index.values, 
         sequence_feature_dataframe(E011)['Total'],
         color=colors)
axs1c.bar(sequence_feature_dataframe(E185B).index.values, 
         sequence_feature_dataframe(E185B)['Total'],
         color=colors)
axs1d.bar(sequence_feature_dataframe(S4_1).index.values, 
         sequence_feature_dataframe(S4_1)['Total'],
         color=colors)
axs1e.bar(sequence_feature_dataframe(S8).index.values, 
         sequence_feature_dataframe(S8)['Total'],
         color=colors)

# individual titles
axs1a.set_title('D133')
axs1b.set_title('E011')
axs1c.set_title('E185B')
axs1d.set_title('S4 1')
axs1e.set_title('S8')

# Set x axis label
fig1.text(0.5, 0.04, 
          'Nucleotide Base', 
          ha='center', 
          va='center',
          fontsize=12)

# Set y axis label
fig1.text(0.06, 0.5, 
          'Amount', 
          ha='center', 
          va='center', 
          rotation='vertical',
          fontsize=12)

# show the plot
plt.show()


# In[22]:


plot_sequence_features(sequence_features_D133)


# In[23]:


length_of_cluster(sequence_features_D133)


# In[24]:


plot_sequence_features(sequence_features_E011)


# In[78]:


dictionary_lengths = {
    'D133' : length_of_code_dataframe(D133)['Species'],
    'E185B' : length_of_code_dataframe(E185B)['Species'],
    'E011' : length_of_code_dataframe(E011)['Species'],
    'S4 1' : length_of_code_dataframe(S4_1)['Species'],
    'S8' : length_of_code_dataframe(S8)['Species']
}

length_df = pd.DataFrame(dictionary_lengths, index=['Nucleotide Base Length'])
length_df


# In[59]:


fig3, axs3 = plt.subplots(1, 1,
                         figsize=(15, 8))

fig3.suptitle('Figure 3 - The length of code for each C.difficle strain')
colors = ['aquamarine', 'mediumaquamarine', 'mediumturquoise',  'lightseagreen', 'teal', 'darkslategrey']
bar_width = 0.4

labels = length_df.index.values
x = np.arange(len(labels)) 

axs3.bar(length_df.index,
         length_df['Species'],
         color=colors,
         width=bar_width)

# x label
axs3.set_xlabel('C. difficile strain',
               fontsize=12)
axs3.set_xticks(labels)

# y label
axs3.set_ylabel('Total Length of Strain (per million)',
               fontsize=12)

axs3.set_yticks(np.arange(0, 6000000, 1000000))

plt.show()


# In[27]:


fig4, (axs4a, axs4b, axs4c, axs4d, axs4e) = plt.subplots(1, 5,
                                                 figsize=(15, 10))
fig4.suptitle('Figure 4 - Percentage of Nucleotide Pairing in each C. difficile strain')

explode = (0.1, 0.1)
colors = ['aquamarine',  'lightseagreen'] 

axs4a.pie(sequence_features_D133[['GC_content', 'AG_content']].T[0], 
          labels=sequence_features_D133[['GC_content', 'AG_content']].T.index.values,
          explode=explode,
          autopct='%1.1f%%',
          colors=colors)

axs4b.pie(sequence_features_E011[['GC_content', 'AG_content']].T[0], 
          labels=sequence_features_E011[['GC_content', 'AG_content']].T.index.values,
          explode=explode,
          autopct='%1.1f%%',
          colors=colors)

axs4c.pie(sequence_features_E185B[['GC_content', 'AG_content']].T[0], 
          labels=sequence_features_E185B[['GC_content', 'AG_content']].T.index.values,
          explode=explode,
          autopct='%1.1f%%',
          colors=colors)

axs4d.pie(sequence_features_S4_1[['GC_content', 'AG_content']].T[0], 
          labels=sequence_features_S4_1[['GC_content', 'AG_content']].T.index.values,
          explode=explode,
          autopct='%1.1f%%',
          colors=colors)

axs4e.pie(sequence_features_S8[['GC_content', 'AG_content']].T[0], 
          labels=sequence_features_S8[['GC_content', 'AG_content']].T.index.values,
          explode=explode,
          autopct='%1.1f%%',
          colors=colors)

# subtitles for each plot
axs4a.set_title('D133')
axs4b.set_title('E011')
axs4c.set_title('E185B')
axs4d.set_title('S4_1')
axs4e.set_title('S8')
 
plt.show()


# In[ ]:




