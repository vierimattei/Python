import csv
import numpy as np
import pandas as  pd
from collections import Counter
from matplotlib import pyplot as plt

data = pd.read_csv('data.csv')
ids = data['Responder_id']
lang_responses = data['LanguagesWorkedWith']

language_counter = Counter()

for response in lang_responses:
    language_counter.update(response.split(';')) 

thelist = language_counter.most_common(15)
        
language = []
popularity = []

for item in thelist:
    language.append(item[0])
    popularity.append(item[1])
  
language.reverse()
popularity.reverse()    

plt.barh(language,popularity)




#plt.style.use('ggplot')
#
#ages_x = np.linspace(25,35,11)
#
#dev_y = np.linspace(38496,73752,11)
#dev_y = dev_y + np.random.random(len(dev_y))*dev_y.min()/10
#
#py_dev_y = np.linspace(45372,83640,11)
#py_dev_y  = py_dev_y + np.random.random(len(py_dev_y ))*py_dev_y .min()/10
#
#ja_dev_y = np.linspace(37810,74583,11)
#ja_dev_y  = ja_dev_y + np.random.random(len(ja_dev_y ))*ja_dev_y .min()/10

#plt.plot(ages_x, py_dev_y, label='Python')
##color='#5a7d9a', marker = '.',linewidth = '3',
#plt.plot(ages_x, ja_dev_y, label='Java')
##color='red', marker = '.', linewidth = '3',
#plt.plot(ages_x, dev_y, label='All')
##color='#444444', linestyle='--', marker = '.', 

#plt.title('salary')
#plt.xlabel('ages')
#plt.ylabel('median salary')
#
#plt.legend()
#
##plt.grid(True) 
#
#plt.tight_layout()
#
#plt.savefig('plot.png')
#
#plt.show()

