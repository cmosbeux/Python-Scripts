#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:18:59 2020

@author: cmosbeux
"""

import textract
import os
from os import listdir
from os.path import isfile, join


#word_file = '../../My_Papers/2020 - Mosbeux - SSH/Ice shelf seasonal flow variations driven by sea surface height.docx'
#word_file = '../../../My_Projects/NASA_Subglacial_Water/Subglacial water dynamics 20:20.docx'
#bibtex = '../../My_Papers/2020 - Mosbeux - SSH/My_Library_June2020.bib'
#biblio_file = 'Biblio-Mosbeux-2020.html'


class cprint:
    @staticmethod
    def green(text):
        return "\033[92m"+text+'\033[0m'
    @staticmethod
    def red(text):
        return "\033[91m"+text+'\033[0m'
    @staticmethod
    def blue(text):
        return "\033[94m"+text+'\033[0m'
        

def display_logo():
    print(cprint.green(" ____ _____ ____ ___ _______ ________   __" ))
    print(cprint.green("|  _ \_   _|  _ \__ \__   __|  ____\ \ / /" ))
    print(cprint.green("| |_) || | | |_) | ) | | |  | |__   \ V / " ))
    print(cprint.green("|  _ < | | |  _ < / /  | |  |  __|   > <  " ))
    print(cprint.green("| |_) || |_| |_) / /_  | |  | |____ / . \ " ))
    print(cprint.green("|____/_____|____/____| |_|  |______/_/ \_\ "))
    print( "\n")
    

def warnings():
    print("\n")
    print(cprint.red("Warnings:"))
    print(cprint.red("---------"))


def extension(file):
    ext = file.split('.')[-1].lower()
    return ext


def select_bibtext(bibpath):
    #look for the bibfiles         
    onlyfiles = [f for f in listdir(bibpath) if isfile(join(bibpath, f))]
    print('\nHere are the bibtex files we found:')
    found = False
    while not found:
        num = 1
        file_list = []
        for file in onlyfiles:
            if extension(file) == 'bib':
                print('\t%d) %s' % (num, file))
                file_list.append(file)
                num+=1
        try:
            file_num = int(input('Select which bibtex file you want to use: '))
        except ValueError:
            print('An integer is required. Please try again.')
        found = True
    #combine path  and file name
    bibpath = bibpath+'/'+file_list[file_num-1] 
    return bibpath
 
    
def select_word(wordpath):
    #look for the bibfiles         
    onlyfiles = [f for f in listdir(wordpath) if isfile(join(wordpath, f))]
    print('\nHere are the bibtex files we found:')
    found = False
    while not found:
        num = 1
        file_list = []
        for file in onlyfiles:
            if extension(file) == 'doc' or extension(file) == 'docx' :
                print('\t%d) %s' % (num, file) )
                file_list.append(file)
                num+=1
        try:
            file_num = int(input('Select which WORD file you want to use: '))
        except ValueError:
            print('An integer is required. Please try again.')
        found = True
    #combine path  and file name
    wordpath = wordpath+'/'+file_list[file_num-1] 
    return wordpath


def get_bibtex():
    print(cprint.green('\nBIBTEX FILE INFO:'))
    print(cprint.green('-----------------'))
    found = False
    while not found:
        try:
            with open('saved_biblio_directory.txt', 'r') as txt:
                saved_biblio = txt.readline()
            if saved_biblio != '':
                print('\nHere is the last bibtex you used : "%s" ' % saved_biblio)
                bibpath = input('If you want to use it press [ENTER] else, enter the path to your BIBTEX file: ')
                if bibpath == '':
                    bibpath = saved_biblio    
                    print('You pressed [Enter]. I will use last bibtex.' )  
                    found = True
                else:
                    try:
                        listdir(bibpath)
                        found = True
                        bibpath = select_bibtext(bibpath)                 
                    except OSError:
                        print('The path seems incorrect... Please enter a new path')               
        except IOError:
            bibpath = input('Enter the path to your BIBTEX file: ') 
            try:
                listdir(bibpath)
                found = True
                bibpath = select_bibtext(bibpath)  
            except OSError:
                print('The path seems incorrect... Please enter a new path')
                
    #save last file:
    with open('saved_biblio_directory.txt', 'w') as txt:
        txt.writelines(bibpath)                           
        
    return bibpath
  
    
def get_file():
    print(cprint.green('\nWORD FILE INFOS:'))
    print(cprint.green('----------------'))
    found = False
    while not found:
        try:
            with open('saved_word_directory.txt', 'r') as txt:
                saved = txt.readline()
            if saved != '':
                print('\nHere is the last word you used : "%s" ' % saved)
                wordpath = input('If you want to use it press [ENTER] else, enter the path to your WORD file: ')
                if wordpath == '':
                    wordpath = saved    
                    print('You pressed [Enter]. I will use last document.') 
                    found = True
                else:
                    try:
                        listdir(wordpath)
                        found = True
                        wordpath = select_word(wordpath)                 
                    except OSError:
                        print('The path seems incorrect... Please enter a new path')                      
        except IOError:
            wordpath = input('Enter the path to your word file: ') 
            try:
                listdir(wordpath)
                found = True
                wordpath = select_word(wordpath)  
            except OSError:
                print('The path seems incorrect... Please enter a new path')
                
    #save last file:
    with open('saved_word_directory.txt', 'w') as txt:
        txt.writelines(wordpath)
    return wordpath



def cite_error(values, valuenames):
    """check if some values are missing for a reference"""
    firstime = True
    for error, errorname in zip(values, valuenames):
        if len(error) == 0:
            if firstime:
                print(cprint.red("(!)  %s (%s)" % (firstauthor,date)))
                firstime = False
            print(cprint.red("\t- missing %s" % (errorname)))


def split(delimiters, string, maxsplit=0):
    import re
    regexPattern = '|'.join(map(re.escape, delimiters))
    return re.split(regexPattern, string, maxsplit)


def title_or_not():
    print('\n')
    a = input('title or not (default y)? [y/n] ') 
    if a == 'n':
        print('\n OK no title... ')
        return False
    else:
        return True
    
def clean_duplicates(line, ld=(' , ', ' : ')):
    for i in ld:
        while ' , ' in line:
            line = line.replace(i, ' ') #remove extra comma from empty data
    return line

class bibliostyle():
    def TC(authors, year, title, editor, booktitle, journal, volume, series, publisher, pages, address, doi):
        """Authors (Last, First. Sec., ...), Title, Journ. Abbrev., Vol, pgs, doi, year"""
        
        bstyle = '{0}: {4}: {2}, {5}, {3}, {6}, {7}, {8}, {1}'
        pstyle = '{0}: {2}, {3}, {4}, {6}, {7}, {1}'

        
        #write a string with the cite in an html format
        #cite_error(info_values, infos)
        if (len(booktitle)>0 or len(publisher)>0) and (len(journal)==0):
            #its a book
            print('book found: %s (%s)' % (authors, year))
            new_line = bstyle.format(authors, year, title, editor, booktitle, series, publisher, address, doi) #trick if address is empty       
        else:
            #its a paper
            new_line = pstyle.format(authors, year, title, journal, volume, number, pages, doi)
            new_line = clean_duplicates(new_line)           
        
        return new_line
    
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

def journal_abrev(journal):
        if journal == 'Geophysical Research Letters':
            journal = 'Geophys. Res. Lett'
        if 'The Cryosphere Discuss' in journal:
            journal = 'Cryosphere Discuss.'
        elif journal == 'The Cryosphere':
            journal = 'Cryosphere'
        if journal == 'Journal of Geophysical Research: Earth Surface':
            journal = 'J. Geophys. Res. Earth Surf.'
        if journal == 'Journal of Geophysical Research: Solid Earth':
            journal = 'J. Geophys. Res. Solid Earth'
        if journal == 'Journal of Glaciology':
            journal = 'J. Glaciol'
        if journal == 'Geoscientific Model Development':
            journal = 'Geosci. Model Dev.'
        if journal == '	Global Planetary Change':
            journal = 'Glob. Planet. Change'
        if journal == 'Proceedings of the National Academy of Sciences':
            journal = 'PNAS'
        if journal == 'Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and Engineering Sciences':
            journal = 'Phil. Trans. R. Soc. A ou Philos. Trans. R. Soc.'
        if journal == 'Nature Geoscience':
            journal = 'Nat. Geosci'
        #journal = '<i>%s</i>' % (journal)
        journal = '%s' % (journal)
        return journal
    
def get_citation_info(ref, papers):
    try:
        firstauthor = papers[ref]['author'].split(',')[0]
    except KeyError:
        print("(!) ref %s:" % (ref))
    authors = papers[ref]['author']
    title = papers[ref]['title']
    
    try:
        date = papers[ref]['year']
    except KeyError:
        print("(!)  %s (?) year missing" % (firstauthor))
        date = ''
    try:                     
        journal = papers[ref]['journal']
        journal = journal_abrev(journal)
    except KeyError:
        journal = ''
    try:
        volume = papers[ref]['volume']
    except KeyError:
        volume = ''
    try:
        number = papers[ref]['number']
        number = '(%s)' % (number)
    except KeyError:
        number = ''
    try:
        pages = papers[ref]['pages']
        pages = pages.replace('--', '&#8212;')
        pages = pages
    except KeyError:
        pages = ''
    try:
        doi =  papers[ref]['doi']
        doi = ', https://doi.org/'+doi
        #doi = ', doi:'+ doi
    except KeyError:
        doi = ''
    #author, title, editor, booktitle, serie, volume, Publisher, place,doi
    try:
        booktitle = papers[ref]['booktitle']
        booktitle = '%s' % (booktitle)
    except KeyError:
        booktitle = ''
    try:
        editor = papers[ref]['editor']
        editor = '%s ' % (editor)
    except KeyError:
        editor = ''
    try:
        series = papers[ref]['series']
        series = '%s' % (series)
    except KeyError:
        series = ''
    try:
        publisher = papers[ref]['publisher']
        publisher = '<i>%s</i>' % (publisher)
    except KeyError:
        publisher = ''
    try:
        address = papers[ref]['address']
        address = '<i>%s</i>' % (address)
    except KeyError:
        address = ''
    
    return title, firstauthor, authors, date, journal, volume, number, pages, doi, booktitle, editor, series, publisher, address
        
#%%    
    
display_logo()

word_file = get_file()  
txt = textract.process(word_file)  
txt = str(txt)  
biblio = ''
bibtex = get_bibtex()
with_title = title_or_not()

with open(bibtex, 'r') as btex:
    for line in btex:
        biblio = biblio + line
        
biblio1 = biblio.split('\n@article')[1:]
biblio2 = []
for i in biblio1:
    bib_temp = i.split('\n@misc')
    for j in bib_temp:
        biblio2.append(j)        
#author, title, editor, booktitle, series, volume, Publisher, place,doi
biblio3 = []
for i in biblio2:
    bib_temp = i.split('\n@book')
    for j in bib_temp:
        biblio3.append(j)
biblio = biblio3        
biblio4 = []
for i in biblio3:
    bib_temp = i.split('\n@inproceedings')
    for j in bib_temp:
        biblio4.append(j)
biblio = biblio4   

biblio5 = []
for i in biblio4:
    bib_temp = i.split('\n@incollection')
    for j in bib_temp:
        biblio5.append(j)
biblio = biblio5

biblio6 = []
for i in biblio5:
    bib_temp = i.split('\n@inbook')
    for j in bib_temp:
        biblio6.append(j)
biblio = biblio6 


#go through the tex file and make a dictionnary of references
papers = {}
infos = ['title', 'author','year', 'journal', 'volume', 'number','pages', 'doi', 'editor', 'booktitle', 'series', 'publisher', 'address']
for item in biblio:
    item_list = item.split('\n')
    paper_ref = item_list[0].strip('{},')
    papers[paper_ref] = {}
    for instance in item_list:
        if 'shorttitle' in instance:
            continue
        for i in infos:
            if i + ' =' in instance:  
                k = instance.split('=')[-1]
                k = k.strip(' {},')
                k = k.replace('{', '')
                k = k.replace('}', '')
                k = k.replace('https://doi.org/', '')
                k = k.replace('doi', '')
                k = k.replace('\xe2\x80\x90', '-')
                #to fix title vs booktitle
                if 'booktitle' in instance:
                    papers[paper_ref]['booktitle'] = k 
                else:
                    papers[paper_ref][i] = k
                    
                    

                
#%%   
             
position = []
new_lines = {}                   
repeats ={}
    
warnings()
papercount, citecount = 0, 0
for ref in sorted(papers.keys()):
    paperfound = False

    #paper info
    title, firstauthor, authors, date, journal, volume, number, pages, doi, booktitle, editor, series, publisher, address = get_citation_info(ref,papers)

    #try to avoid crash but maybe not necessary...
    try:
        #get authors initial only
        authors = authors.split(' and ')
        authors_new = []
        x=20
        for author in authors:
                lastname = author.split(',')[0]
                try:
                    firstnames = author.split(',')[1].split(' ')[1:]
                except IndexError:
                    print("(!)  %s (%s) does not have first name..." % (firstauthor,date))
                for i in range(len(firstnames)):
                    try:
                        firstnames[i] = firstnames[i][0] + '.'
                    except IndexError:
                        print("(!)  %s (%s) problem with first name" % (firstauthor,date))
                author_new = lastname + ', ' + ' '.join(firstnames)
               
                if len(authors) > 2:
                    authors_new.append(author_new)
                    authors = ' '+', '.join(authors_new[:-1]) + ' and ' + authors_new[-1]
            
                if len(authors_new) > x:
                	authors = ', '.join(authors_new[:x-1]) + ' et al.' 
                elif len(authors) == 2:
                    authors_new.append(author_new)
                    authors = authors_new[0] + ' and ' + authors_new[-1]
                elif len(authors) == 1:
                    authors = author_new
    except KeyError:
        'dummy'

    info_values = [title, authors, date, journal, volume, number, pages, doi]
    
    #deal with special characters
    firstauthor_unicode = firstauthor
    if 'ü' in firstauthor:
        firstauthor_unicode = firstauthor_unicode.replace('ü', '\\xc3\\xbc')        
        
    if firstauthor_unicode in txt:
        year = 0000
        #indexes = [n for n in range(len(txt)) if txt.find(firstauthor, n) == n] 
        indexes = list(find_all(txt, firstauthor_unicode))

        for index in indexes:
            for i in range(35): #limit of 40 caracters after the name to look after to get the date
                yearscan = txt[index+i:index+i+4]
                if yearscan == date:
                    #paperfound = True
                    break
            if yearscan == year: 
                citecount+=1
            elif yearscan.isdigit():
                paperfound = True
                papercount+=1
                year = yearscan
                position = index #save only the first time the reference was cited
    
    if paperfound:
        #write a string with the cite in an html format
        cite_error(info_values, infos)
        #make new line entry with correct bibliostyle
        new_line = bibliostyle.TC(authors, year, title, editor, booktitle, journal, volume, series, publisher, pages, address, doi)
        new_lines[position] = new_line
        

print(cprint.green('\n %d cites found for %d unique references.' % ( citecount, papercount)))

#%write the bibliography dictionnary to an html document
with open('Bibliography.html', 'w') as f:   
    k=1
    alphabetical = input('Alphabetical or text order ? [1:default / 2] ')     
    
    #start writing the html output file
    #f.write('</p>')
    if alphabetical.lower() == '2':
       num_sort = list(new_lines.keys())
       num_sort.sort()
       for pos in num_sort: 
           f.write('<p>')
           full_line = '<b>  ['+ str(k)+ ']</b> ' +new_lines[pos]
           f.writelines(full_line)
           f.write('</p>')
           k+=1
    else:
       num_sort = new_lines.values()
       #num_sort.sort()
       for pos in num_sort:
           f.write('<p>&emsp;')
           full_line = pos 
           f.writelines(full_line)
           f.write('.</p>')
           f.write('\n')
           k+=1

    #f.write('</p>')

print('\nBibliography have been saved in: ')
os.system("pwd")
open_url = input ("Open URL ?  [y/n] ")
if open_url.lower() == 'y':
    os.system("\nopen Bibliography.html")
           
#%%
# #its a book
#     if len(address)>0:
#         if with_title:
#             new_line = '{0} ({1}) {2}. {3}{4}, {5}, {6}: {7}{8}'.format(authors, year, title, editor, booktitle, series, publisher, address, doi)
#         else:
#             new_line = '{0} ({1}) {2}. {3}{4}, {5}, {6}: {7}{8}'.format(authors, year, title, editor, booktitle, series, publisher, address, doi)
#     else:
#         if with_title:
#             new_line = '{0} ({1}) {2}. {3}{4}, {5}, {6}{7}'.format(authors, year, title ,editor, booktitle, series, publisher, doi)
#         else:
#             new_line = '{0} ({1}) {2}{3}, {4}, {5}{6}'.format(authors, year, editor, booktitle, series, publisher, doi)
#     print('book found: %s (%s)' % (authors, year))
# else:
# #its a paper
        
#     if len(pages)>0:
#         if with_title:
#             new_line = '{0} ({1}) {2}. {3}, {4}{5}, {6}{7}'.format(authors, year, title, journal, volume, number, pages, doi)
#         else:
#             new_line = '{0} ({1}) {2}{3}  '.format(authors, year, journal, doi)
#     else:
#         if with_title:
#             new_line = '{0} ({1}) {2}. {3}{4}{5}, {6}'.format(authors, year, title, journal, volume, number, doi)
#         else:
#             new_line = '{0} ({1}) {2}{3}  '.format(authors, year, journal, doi)
