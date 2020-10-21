#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np

#Increasing the width of the notebook (visual difference only)
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

import pandas as pd

from MONDquantities import*

#All the measurements use a reduced Hubble constant. As the Hubble constant can be expressed in terms
#of Miglrom's constant and speed of light, we express it as that. I modified a0 = 1.14 * 10**-10 so 
#that H0 is now very close to the latest (upper) estimate. Now a0 and H0 are consistent.
H0 = 2*pi*(a0)/c

#The Hubble constant used is h50, given as H0/50 km/s/Mpc. We need to multiply most of the quantities in
#the table by this to obtain 
h50 = H0/(50*10**3/(1000*kp))

#In the database, there is an optional parameter c that indicates a different way of calculating the
#error. Only some lines have it, so it throws pandas off. Removed it in the file itself
#a2065, a2063, ngc5846 are missing from reiprich 2001 for rho_0. Had to add them in by hand 
#in the rho_0 file.
#In the PhD thesis Reiprich cut the table wrong so they're invisible! Unbelievable

#Importing the text file
df = pd.read_csv('cluster_database/cluster_data.txt', delimiter='\s+', header = None)

#As the columns have no names in the initial file, we add the names here
df.columns = ['name', 'beta_frame', 'beta+', 'beta-', 'r_c_frame', 'r_c+', 'r_c-', 'T', 'T+', 'T-',
              'm5_frame', 'm5+', 'm5-', 'r5', 'r5+', 'r5-', 'm2', 'm2+', 'm2-', 'r2', 'r2+', 'r2-',
              'mtot_frame', 'ref']

# Using readline to open file containing rho_0 values and storing each line in Lines
file1 = open('cluster_database/rho_0.txt', 'r') 
Lines = file1.readlines() 
file1.close()

#Initialising numpy array to hold the values of rho_0
rho_0_np = np.zeros((len(Lines),1))

#Storing the content of each line into the numpy array
for i, line in enumerate(Lines):
    
    #assigning each line to a member of the rho_0 numpy array
    rho_0_np[i] = line

#Flattening the array before putting it in the dataframe    
rho_0_np = rho_0_np[:,0]
#Adding the rho_0 values to the dataframe
df['rho_0_frame'] = pd.Series(rho_0_np, index = df.index)
    
#The radii are given in units of kpc/h50. To get them in kpc we need to divide by h50
#Only one useful for now to be scaled is r_c
df.loc[:, 'r_c_frame'] = df.loc[:, 'r_c_frame']/h50

# Displaying the full, uncut database with all parameters
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(df.head(5))


# In[2]:


h50

