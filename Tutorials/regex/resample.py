# -*- coding: utf-8 -*-
import re

text_to_search = '''
abcdefghijklmnopqrstuvwxyz
ABCDEFGHIJKLMNOPQRSTUVWXYZ
1234567890

Ha HaHa

coreyms.com

MetaCharacters (Need to be escaped)
. ^ $ * + ? { } [ ] \ | ( )

321-555-4321
123.555.1234
123*555*1234
800-555-1234
900-555-1234

Mr. Schafer
Mr Smith
Ms Davis
Mrs. Robinson
Mr. T

cat
mat
pat
bat

'''

sentence = 'Start a sentence and then bring it to an end'

#Compiling string into a regex
pattern = re.compile(r'M(r|rs)\.?\s[A-Z]\w*')

matches = pattern.finditer(text_to_search)

for match in matches:
    print(match)
    
# print(text_to_search[1:4])

# with open('faked.txt','r') as faked:
#     contents = faked.read()
    
#     matches = pattern.finditer(contents)
    
#     for match in matches:
#         print(match)