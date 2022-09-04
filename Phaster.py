#!/usr/bin/env python
# coding: utf-8

# # Libraries

# In[1]:


import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt


# In[2]:


import requests
import wget
# importing webbrowser python module
import webbrowser


# In[3]:


genome_summary = pd.read_excel('Genome summary for SJByrne.xlsx', header=3)
display(genome_summary)


# # PHASTER

# In[5]:


#Assigning URL to be opened
URL = "https://phaster.ca/"

#Open url in default browser
webbrowser.open(URL, 
                new=2) # new = 2, open URL in same tab.


# In[6]:


def import_fasta_files(file):
    '''
    
    '''
    # create empty list to paste the fasta data into
    fasta_data = []
    with open(file, 'r') as fasta_file:
        fasta_data = fasta_file.read()
        fasta_data = fasta_data.strip().split("\n")

    return fasta_data


# In[7]:


def API_checklist(api):
    '''
    This function:
    
    Parameters:
    -----------
    '''
    result = api.reason
    print(f'result = {result}')
    
    text = api.text
    print(f'text = {text}')
    
    return


# In[8]:


def phaster_api(data):
    
    # 1. connect with Phaster API
    response_API = requests.post(url="https://phaster.ca/phaster_api",
                                json=data)

    # 2. checking status code of the API
    if response_API.status_code == 200:
        print("Succesful connection with API.")
    elif response_API.status_code == 404:
        print("Unable to reach URL.")
    else:
        print("Unable to connect API or retrieve data.")
    
    # 3. run the previous function
    print(API_checklist(response_API))
    
    # 4. Find accession number
    accession_number = response_API.text.replace('"', ' ').split()[3]
    
    return accession_number


# In[9]:


import json 

data = import_fasta_files('S8_genome.fasta')


# ---------------------

# In[10]:


response_API = requests.post(url="https://phaster.ca/phaster_api")

# checking status code of the API
if response_API.status_code == 200:
    print("Succesful connection with API.")
elif response_API.status_code == 404:
    print("Unable to reach URL.")
else:
    print("Unable to connect API or retrieve data.")


# 200 : OK. It means we have a healthy connection with the API on web.
# 
# GET: retrieve information (like search results). This is the most common type of request. Using it, we can get the data we are interested in from those that the API is ready to share.
# 
# POST: adds new data to the server. Using this type of request, you can, for example, add a new item to your inventory.
# 
# An API (Application Programming Interface) is a set of rules that are shared by a particular service. These rules determine in which format and with which command set your application can access the service, as well as what data this service can return in the response. The API acts as a layer between your application and external service. You do not need to know the internal structure and features of the service, you just send a certain simple command and receive data in a predetermined format.

# In[11]:


API_checklist(response_API)


# In[12]:


# find accession number
accession_number = response_API.text.replace('"', ' ').split()[3]
accession_number


# ---------------------

# # PHASTER Results

# In[4]:


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


# ### Functions

# In[44]:


def phaster_results(accession_number):
    
    # 1. 
    response_phasterAPI = requests.get(url=f'https://phaster.ca/phaster_api?acc={accession_number}')
    
    # 2.
    url = f'http://phaster.ca/submissions/{accession_number}.zip'

    # Downloading the file by sending the request to the URL
    req = requests.get(url)
 
    # Split URL to get the file name
    filename = url.split('/')[-1]
 
    # Writing the file to the local file system
    with open(filename,'wb') as output_file:
        output_file.write(req.content)
    
    # specifying the zip file name
    file_name = f'{accession_number}.zip'
  
    # opening the zip file in READ mode
    with ZipFile(file_name, 'r') as zip:
    # printing all the contents of the zip file
        zip.printdir()
    # extracting all the files
        zip.extractall() 
        
    return 


# I renamed the files (details, summary, phaster.fna) with the genome name in front of the file externally

# In[56]:


phaster_results(phaster_df['Accession Number'][0])
print(phaster_df['Genome'][0])


# In[57]:


phaster_results(phaster_df['Accession Number'][1])
print(phaster_df['Genome'][1])


# In[58]:


phaster_results(phaster_df['Accession Number'][2])
print(phaster_df['Genome'][2])


# In[59]:


phaster_results(phaster_df['Accession Number'][3])
print(phaster_df['Genome'][3])


# In[60]:


phaster_results(phaster_df['Accession Number'][4])
print(phaster_df['Genome'][4])


# ------

# # Summary Files

# ### Pre-processing Functions

# In[63]:


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


# In[68]:


def summaryname_1(summaryname):
    '''
    
    '''
    
    summary = open(summaryname).read().split()
    
    del summary[0:329]
    
    summary.pop(17)
    
    summary_df = summary_to_dataframe(summary)
    
    return summary_df


# In[85]:


def summaryname_2(summaryname):
    '''
    
    '''
    summary = open(summaryname).read().split()
    
    del summary[0:336]
    summary.pop(17)
    
    summary_df = summary_to_dataframe(summary)
    return summary_df


# ### Genome Summary Files

# In[69]:


S8_summary = summaryname_1('S8_Genome_summary.txt')
display(S8_summary)


# In[70]:


S4_1_summary = summaryname_1('S8_Genome_summary.txt')
display(S4_1_summary)


# In[86]:


E185B_summary = summaryname_2('E185B_Genome_summary.txt')
display(E185B_summary)


# In[87]:


E011_summary = summaryname_2('E011_Genome_summary.txt')
display(E011_summary)


# In[88]:


D133_summary = summaryname_2('D133_Genome_summary.txt')
display(D133_summary)


# ------

# # BLAST Query Pre-processing

# ### Functions

# In[247]:


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


# In[323]:


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
           
   


# ### Genome Files

# In[242]:


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

# In[248]:


D133_nucleotide_position = nucleotide_position(D133_summary['REGION_POSITION'])
D133_nucleotide_position


# In[ ]:


D133_phages = find_phages(D133_nucleotide_position, D133)


# #### E011

# In[239]:


E011_nucleotide_position = nucleotide_position(E011_summary['REGION_POSITION'])
E011_nucleotide_position


# In[1]:


E011_phages = find_phages(E011_nucleotide_position, E011)


# #### E185B

# In[240]:


E185B_nucleotide_position = nucleotide_position(E185B_summary['REGION_POSITION'])
E185B_nucleotide_position


# In[ ]:


E185B_phages = find_phages(E185B_nucleotide_position, E185B)


# #### S4 1

# In[241]:


S4_1_nucleotide_position = nucleotide_position(S4_1_summary['REGION_POSITION'])
S4_1_nucleotide_position


# In[ ]:


S4_1_phages = find_phages(S4_1_nucleotide_position, S4_1)


# #### S8

# In[238]:


S8_nucleotide_position = nucleotide_position(S8_summary['REGION_POSITION'])
S8_nucleotide_position


# In[249]:


S8 = import_fasta_files('S8_genome.fasta')
display(S8)


# In[331]:


S8_phages = find_phages(S8_nucleotide_position, S8)
display(S8_phages)


# ---

# # BLAST

# ### Specific Libraries

# In[116]:


import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import GenBank 
from Bio.SeqIO.FastaIO import SimpleFastaParser

# to use BLAST library
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# ## BLAST over the Internet

# In[ ]:


help(NCBIWWW.qblast) 


# In[ ]:


def blast_over_the_internet(sequence):
    result_handle = NCBIWWW.qblast(
                    program="blastn",               # the program to search
                    database="nt",                  # Which database to search against 
                    sequence=sequence,              # the sequence to search
                    
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


# demonstrating doing a function for BLAST takes forever:
# 
# Running BLAST locally (as opposed to over the internet, see Section ‍7.1) has at least major two advantages:
# 
# Local BLAST may be faster than BLAST over the internet;
# Local BLAST allows you to make your own database to search for sequences against.
# Dealing with proprietary or unpublished sequence data can be another reason to run BLAST locally. You may not be allowed to redistribute the sequences, so submitting them to the NCBI as a BLAST query would not be an option.
# 
# Unfortunately, there are some major drawbacks too – installing all the bits and getting it setup right takes some effort:
# 
# Local BLAST requires command line tools to be installed.
# Local BLAST requires (large) BLAST databases to be setup (and potentially kept up to date).
# To further confuse matters there are several different BLAST packages available, and there are also other tools which can produce imitation BLAST output files, such as BLAT.

# In[ ]:


blast_over_the_internet(sequence)


# In[ ]:


get_ipython().run_line_magic('timeit', 'blast_over_the_internet(sequence)')


# We use the function qblast() in the Bio.Blast.NCBIWWW module to call the online version of BLAST. This has three non-optional arguments:
# 
# The first argument is the blast program to use for the search, as a lower case string. The options and descriptions of the programs are available at https://blast.ncbi.nlm.nih.gov/Blast.cgi. Currently qblast only works with blastn, blastp, blastx, tblast and tblastx.
# The second argument specifies the databases to search against. Again, the options for this are available on the NCBI Guide to BLAST ftp://ftp.ncbi.nlm.nih.gov/pub/factsheets/HowTo_BLASTGuide.pdf.
# The third argument is a string containing your query sequence. This can either be the sequence itself, the sequence in fasta format, or an identifier like a GI number.
# The qblast function also take a number of other option arguments which are basically analogous to the different parameters you can set on the BLAST web page. We’ll just highlight a few of them here:
# 
# The argument url_base sets the base URL for running BLAST over the internet. By default it connects to the NCBI, but one can use this to connect to an instance of NCBI BLAST running in the cloud. Please refer to the documentation for the qblast function for further details.
# The qblast function can return the BLAST results in various formats, which you can choose with the optional format_type keyword: "HTML", "Text", "ASN.1", or "XML". The default is "XML", as that is the format expected by the parser, described in section ‍7.3 below.
# The argument expect sets the expectation or e-value threshold.

# ## BLAST Locally

# From within Biopython we can use the NCBI BLASTX wrapper from the Bio.Blast.Applications module to build the command line string, and run it:

# The BLAST+ package for running BLAST locally was downloaded at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# 
# using the source link and code --> ncbi-blast-2.13.0+-win64.exe 

# In[327]:


from Bio.Blast.Applications import NcbiblastxCommandline


# In[328]:


help(NcbiblastxCommandline)


# In[ ]:


# create a blastx query to compare nucleotide to nucleotide database 
# set the e value cutoff at 1e-20 
# choose xml output format 
# limit the output to the best BLAST match for each query sequence 


# In[361]:


# what i would put in the command line
blastx –query fastafilename.fasta –db database to query -evalue 1e-20 –outfmt 5 -num_descriptions 1 -num_alignments 1


# In[351]:


blastx_cline = NcbiblastxCommandline(query=S8_phages['phage_code'][0],                            #          
                                     db="",                                  # the downloaded database
                                     evalue=0.001,
                                     outfmt=5, 
                                     out="opuntia.xml")


# In[356]:


stdout, stderr = blastx_cline()


# In[ ]:


E_VALUE_THRESH = 1e-20
for record in NCBIXML.parse(open("my_blast.xml")):
    if record.alignments : #skip queries with no matches
        print "QUERY: %s" % record.query[:60]
    for align in record.alignments:
        for hsp in align.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print "MATCH: %s " % align.title[:60]
                print hsp.expect


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ---

# ------

# # Phage Region Files

# ### Functions

# In[26]:


def character_present(character, test_list):
    '''
    
    '''
    for element in character:
        if element in test_list:
            return 0
    return 1


# In[27]:


def phage_region_to_dataframe(file):
    '''
    
    '''
    
    # 1. pieces is the total length of the list/2
    pieces = len(file)/2
    
    # 2. split the list into format of [consensus number, code] 
    fasta_file_split = np.array_split(file, pieces)
    
    # 3. make a dataframe
    summary_df = pd.DataFrame(fasta_file_split,
                             columns=['Phage Region', 'Code']) 
    
    return summary_df


# In[52]:


def phage_region(filename):
    
    phage_region = open(filename).read().split()
    
    character = '>'

    phage_region_result = [element for element in phage_region if character_present(element, character)]
    
    return phage_region_result


# ### Genome Phaster Region Files

# In[325]:


S8_phage_region = phage_region('S8_Genome_phage_regions.fna')


# In[110]:


S4_1_phage_region = phage_region('S4_1_Genome_phage_regions.fna')


# In[ ]:


E185B_phage_region = phage_region('E185B_Genome_phage_regions.fna')


# In[ ]:


E011_phage_region = phage_region('E011_Genome_phage_regions.fna')


# In[ ]:


D133_phage_region = phage_region('D133_Genome_phage_regions.fna')


# ## Detail Files

# ### Functions

# In[ ]:


def character_present(character, test_list):
    '''
    
    '''
    for element in character:
        if element in test_list:
            return 0
    return 1


# ### Genome Files

# In[113]:


detail = open("S8_Genome_detail.txt").read()#.strip().split()
detail


# In[ ]:


S8_detail = open("S8_Genome_detail.txt").read().strip().split()

S4_1_detail = open("S8_Genome_detail.txt").read().strip().split()

S8_detail = open("S8_Genome_detail.txt").read().strip().split()

S8_detail = open("S8_Genome_detail.txt").read().strip().split()

S8_detail = open("S8_Genome_detail.txt").read().strip().split()


# ------

# In[ ]:


# close connection with the server
response_API.close()
------
response_phasterAPI.close()


# ------

# # MACHINE LEARNING

# ### Libraries

# In[373]:


from Bio import SeqIO
from Bio import Entrez


# In[389]:


# Install the tensorflow addons package,
get_ipython().system('pip install tensorflow-addons')


# In[390]:


# downloading the plugin profile because it wasnt working properly 
get_ipython().system('pip install -U tensorboard_plugin_profile')


# In[400]:


get_ipython().system('pip3 install --upgrade  tensorflow')


# In[407]:


# Import modules
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential

from keras.callbacks import EarlyStopping
from keras.callbacks import ReduceLROnPlateau

from datetime import datetime

import tensorflow_io as tfio


# In[403]:


# printing the version of tensorflow
print("TensorFlow version: ", tf.__version__)


# In[405]:


# confirming that tensorflow can access the GPU

device_name = tf.test.gpu_device_name()
if not device_name:
    raise SystemError('GPU device not found')
print('Found GPU at: {}'.format(device_name))


# ---

# In[388]:


# Importing the Tn6215 genome into python

# Have to set my email otherwise Entrez does not work
Entrez.email = 'sjjbyrne@hotmail.com'

handle = Entrez.efetch(db="nuccore", id="KC166248", rettype="gb", retmode="text")
genome = SeqIO.parse(handle, "genbank")
for record in genome:
    print(f'     The record id = \n {record.id} \n \n     The record descriptions = \n {record.description} \n \n     The record = \n {record} \n \n     The record sequence = \n {record.seq}')


# ---

# ### 1. One-Hot Encoder

# In[ ]:


# one-hot encoder on the Tn6215 genome

one_hot = tfio.genome.sequences_to_onehot(record.seq)
print(one_hot)
print(one_hot.shape)


# In[ ]:


# one hot encoder for genome files


# In[ ]:


def create_model():
    model = Sequential([
        layers.Input(),
        layers.Conv2D(),
        layers.MaxPooling2D(),
        layers.Flatten(),
        layers.Dense()
        name='CNN_model'
    ])
    
    return model


# In[ ]:


model_1 = create_model()

model_1.summary()


# In[ ]:


# Create a TensorBoard callback
logs = "logs/" + datetime.now().strftime("%Y%m%d-%H%M%S")

tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=logs,
                                                      histogram_freq=1,
                                                      profile_batch='21, 31')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ## DNA Sequence Similarity

# #### The Needleman-Wunsch Algorithm

# In[364]:


def nw(x, y, match = 1, mismatch = 1, gap = 1):
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


# In[ ]:


x = 
y = 


# In[363]:


print(nw(x, y, gap = 1))

# gap = 


# In[ ]:


# timeit 


# In[ ]:


# accuracy compared to BLAST results


# #### 

# In[ ]:





# In[ ]:





# ## Classification Algorithms

# Classification is one of the most studied tasks in machine learning. The principle of classification is based on the predicted attribute to predict the class of the target attribute specified by the user. In genomics, the key issues are genome classification and sequence annotation. In the mining of biological sequences, widely used algorithms include:
# - fuzzy sets, 
# - neural networks, 
# - genetic algorithms, 
# - rough sets. 
# 
# There are also many general classification models, such as:
# - naive Bayesian networks, 
# - decision trees, 
# - neural networks,
# - rule learning using evolutionary algorithms.

# In[ ]:


x = 
y = 


# In[ ]:


from sklearn.model_selection import train_test_split

x_train, x_test, y_train, y_test = train_test_split(x, y , 
                                                    test_size=0.25, 
                                                    random_state=0)


# In[ ]:


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


# In[ ]:


from sklearn.metrics import accuracy_score, precision_score, recall_score

accuracy, precision, recall = {}, {}, {}

for key in models.keys():
    
    # Fit the classifier model
    models[key].fit(X_train, y_train)
    
    # Prediction 
    predictions = models[key].predict(X_test)
    
    # Calculate Accuracy, Precision and Recall Metrics
    accuracy[key] = accuracy_score(predictions, y_test)
    precision[key] = precision_score(predictions, y_test)
    recall[key] = recall_score(predictions, y_test)


# In[ ]:


df_model = pd.DataFrame(index=models.keys(), columns=['Accuracy', 'Precision', 'Recall'])
df_model['Accuracy'] = accuracy.values()
df_model['Precision'] = precision.values()
df_model['Recall'] = recall.values()

df_model


# In[ ]:


from sklearn.metrics import confusion_matrix

cm = confusion_matrix(y_test, predictions)

TN, FP, FN, TP = confusion_matrix(y_test, predictions).ravel()

print('True Positive(TP)  = ', TP)
print('False Positive(FP) = ', FP)
print('True Negative(TN)  = ', TN)
print('False Negative(FN) = ', FN)

accuracy =  (TP+TN) /(TP+FP+TN+FN)

print('Accuracy of the binary classification = {:0.3f}'.format(accuracy))


# In[ ]:


ax  = df_model.plot.bar(rot=45)
ax.legend(ncol= len(models.keys()), 
          bbox_to_anchor=(0, 1), 
          loc='lower left', 
          prop={'size': 14})
plt.tight_layout()


# In[ ]:





# ## Clustering Algorithms

# In[ ]:





# ## something else

# In[366]:


import PyMix


# In[367]:


pip install PyMix


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ------
