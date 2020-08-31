#Taking the information on the Virgo Cluster galaxies. The data comes from 
#http://www.atlasoftheuniverse.com/galgrps/vir.html

#Importing pandas for data acquisition/manipulation
import pandas as pd

#Importing re to use regular expressions
import re

#defining an empty list galaxy[] to contain each line of our file
galaxies = []

#Using a contextmanager through 'with' so after we are done with the file it 
#gets closed automatically
#Copying the data from the original file so we can modify that and leave the 
#original data intact. Opening 'data.csv' file equates to creating it if it
#doesnt yet exist. Here we'' copy the formatted data from galaxies

with open('data.txt', 'r') as data:
    with open('data.csv', 'r+') as data_csv:
        
        #writing data from the original file to the galaxy list
        #line by line so we don't use up all the memory at once
        for line in data:
            galaxies.append(line)
                        
#Contexts for data and data_csv are now closed. 
            
# defining the re patterns we want to find. r means include all characters,
# e.g. a newline is coded as \n. 

#The first few lines have a space before the newline, so that messes up the
#next regex. So we look for space+newline so we can get rid of the space so that
#all lines have the same structure
pattern_1 = re.compile(r' \n')

#Finding any number of at least two spaces
#so we avoid considering the newlines! (there is always only one newline, 
#more than one space, only one space in the same field (e.g. for coordinates))
pattern_2 = re.compile(r'\s\s+')

#Defining a new list galaxies_csv for the csv format
galaxies_csv = []

i = 0
#finding the matches in the galaxies
for row in galaxies:
    
    #Fixing the space+newline in the first few lines
    row = re.sub(pattern_1, '\n', row)
    #Directly subbin ',' for each match to the regex pattern in each row
    row = re.sub(pattern_2, ',', row)
    #Overwriting the value of each row with the subbed expression
    galaxies_csv.append(row)
    #Incrementing counter after subbin so we start from 0
    i +=1

#Notifying whether lists have same length (which has to be the case for them to
#be correct))
if len(galaxies) == len(galaxies_csv):
    print('Lists have same length, comversion should be OK')
else:
    print('Lists have different length, error in conversion')

with open('data.csv', 'w') as data_csv:
    for row in galaxies_csv:
        data_csv.writelines(row)




#We open another context to 
#manipulate the copy we created to make it into a csv file. Using r+ which
#means we can read and write to it








    
