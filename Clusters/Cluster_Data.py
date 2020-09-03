#Taking the information on the Virgo Cluster galaxies. The data comes from 
#http://www.atlasoftheuniverse.com/galgrps/vir.html

#importing numpy to work with arrays etc
import numpy as np

#Importing pandas for data acquisition/manipulation
import pandas as pd

#Importing re to use regular expressions
import re

#Importing astropy to use constants and coordinates, units to do conversions,
#WMAP cause we need the Hubble constant, SkyCoord to convert to cartesian
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo
from astropy.coordinates import SkyCoord

#importing functions I made that are in this directory
import my_functions as myfun

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

#the literal newline characters at the end of each line
pattern_3 = re.compile(r'\n')

#matching a single comma at the beginning of a line
pattern_4 = re.compile(r'^,')

#Defining a new list galaxies_csv for the csv format
galaxies_csv = []

i = 0
#finding the matches in the galaxies
for row in galaxies:
    
    #Fixing the space+newline in the first few lines
    row = re.sub(pattern_1, '\n', row)
    #Directly subbin ',' for each match to the regex pattern in each row
    row = re.sub(pattern_2, ',', row)
    #removing the literal \n for newline at the end of each row
    row = re.sub(pattern_3, '', row)
    #removing the comma at the beginning of lines
    row = re.sub(pattern_4, '', row)
    #Overwriting the value of each row with the subbed expression
    galaxies_csv.append(row)
    #Incrementing counter after subbin so we start from 0
    i +=1

#Notifying whether lists have same length (which has to be the case for them to
#be correct))
if len(galaxies) == len(galaxies_csv):
    print('Lists have same length, conversion should be OK')
else:
    print('Lists have different length, error in conversion')

#Removing the 0th, 2nd and 3rd elements of the list, which are redundant. After
#deleting [0], elements shift so we have to delete elements 1 and 2 instead
del galaxies_csv[0]
del galaxies_csv[1:3]

#Re-formatting the 1st line with correct names for each of the 0 fields, with
#unit names matching those in astropy
galaxies_csv[0] = 'Name,Right_Ascension,Declination,Magnitude,Type,Size(angle),\
Size(klyr),Recession_Velocity(km/s),Other_Names'

with open('data.csv', 'w') as data_csv:
    for row in galaxies_csv:
        data_csv.write(row + '\n')

#reading in the formatted csv file we made as a data frame
galaxies_frame = pd.read_csv('data.csv')

#Removing redundant series (columns) from the data frame: the size in arcminutes,
#the Other_Names and the type field and putting into a new frame
galaxies_info = galaxies_frame.drop(['Size(angle)','Other_Names','Type'], axis=1)

#Changing the unit of the dataframe using my functions with astropy, then 
#renaming the columns of the data frame to reflect the change in unit
galaxies_info['Size(klyr)'], size_unit = myfun.convert_unit(galaxies_frame['Size(klyr)'],u.klyr,u.kpc)
myfun.update_series_name(galaxies_info,'Size(klyr)', size_unit)

#converting recession velocity to distance using astropy
galaxies_info['Recession_Velocity(km/s)'], dist_unit = myfun.convert_vel_dist(galaxies_info['Recession_Velocity(km/s)'],u.km, u.s)
dist_name = 'Distance'+'('+dist_unit +')'
galaxies_info.rename(columns={'Recession_Velocity(km/s)':dist_name}, inplace=True)


