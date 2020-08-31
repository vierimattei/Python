#Taking the information on the Virgo Cluster galaxies. The data comes from 
#http://www.atlasoftheuniverse.com/galgrps/vir.html

#Importing pandas for data acquisition/manipulation
import pandas as pd

#Importing re to use regular expressions
import re

#Using a contextmanager through 'with' so after we are done with the file it 
#gets closed automatically

#Copying the data from the original file so we can modify that and leave the 
#original data intact

with open('data.txt', 'r') as data:
    with open('data.csv', 'w') as data_csv:
        
        for line in data:
            data_csv.write(line)

#Contexts for data and data_csv are now closed. We open another context to 
#manipulate the copy we created to make it into a csv file. Using r+ which
#means we can read and write to it

with open('data.csv', 'r+') as data_csv:
    galaxies = []
    
    #galaxies is initially an empty list, then every line from the text file
    #is appended to it, so e.g. galaxy[0] is the 1st line of the file etc.
    for line in data_csv:
        galaxies.append(line)

# defining the re patter we want to find. r means include all characters,
# e.g. a newline is coded as \n. Finding any number of at least two spaces
#so we avoid considering the newlines!
pattern = re.compile(r'\s\s+')

i = 0
#finding the matches in the galaxies
for row in galaxies:
    matches = pattern.finditer(row)
    i += 1
    print('Row ' + str(i))
    j = 0
    # matches.append(pattern.finditer(row))
    for match in matches:
        j += 1
        match_pos = (match.span())
        matched_str = row[match.span()[0]:match.span()[1]]
        print(matched_str)
        print('Match ' + str(j) + '=' + '"' + matched_str + '"')
        row.replace(matched_str, ',')
        print(row)


# for element in galaxies:
#     matches.append(pattern.finditer(galaxies[0]))



# for match in matches:
#     print(match)
#     if match != '\n':
#         print('not a newline!')
#     else:
#         print('newline')

    
