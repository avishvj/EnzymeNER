# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 10:35:15 2021

@author: one
"""
#import numpy as np
#import pandas as pd
import json
import re
import os
import time
#from fuzzywuzzy import fuzz

'''
file_path = 'D:/BMR-DS/Project 1/Processing/PMC7093903_bioc.json'
with open(file_path, 'r', encoding = 'utf-8') as file: # windows 10 has the 'gbk' codec problem without encoding = 'utf-8'
    fulltext = json.load(file)                         # fulltext is a dict in python, after using json.load().
file.close()
'''

def enzyme_ase(word):
    matchword = ''
    first_match = ''
    enzymes_list = []
    flag = 1 # If flag==0, search_list only has first_match.
    if re.match('(.*)synthase', word):
        first_match = 'synthase'
        flag = 0
    elif re.match('(.*)lyase', word):
        first_match = 'lyase'
        flag = 0
    elif re.match('(.*)ligase', word):
        first_match = 'ligase'
        flag = 0
    elif re.match('(.*)ease', word):
        first_match = 'ease'
        if re.match('(.*)permease', word):
            matchword = 'permease'
        elif re.match('(.*)protease', word):
            matchword = 'protease'
        else:
            matchword = 'other'
    elif re.match('(.*)lase', word):
        first_match = 'lase'
        if re.match('(.*)methylase', word):
            matchword = 'methylase'
        elif re.match('(.*)horylase', word):
            matchword = 'horylase'
        elif re.match('(.*)cyclase', word):
            matchword = 'cyclase'
        elif re.match('(.*)hydrolase', word):
            matchword = 'hydrolase'
        elif re.match('(.*)oxylase', word):
            matchword = 'oxylase'
        else:
            matchword = 'other'
    elif re.match('(.*)rase', word):
        first_match = 'rase'
        if re.match('(.*)transferase', word):
            matchword = 'transferase'
        elif re.match('(.*)saturase', word):
            matchword = 'saturase'
        else:
            matchword = 'other'
    elif re.match('(.*)tase', word):
        first_match = 'tase'
        if re.match('(.*)synthetase', word):
            matchword = 'synthetase'
        elif re.match('(.*)reductase(.*)', word):
            matchword = 'reductase'
        else:
            matchword = 'other'
    elif re.match('(.*)idase', word):
        first_match = 'idase'
        if re.match('(.*)oxidase', word):
            matchword = 'oxidase'
        elif re.match('(.*)peptidase', word):
            matchword = 'peptidase'
        else:
            matchword = 'other'
    elif re.match('(.*)nase', word):
        first_match = 'nase'
        if re.match('(.*)proteinase', word):
            matchword = 'proteinase'
        elif re.match('(.*)oxygenase', word):
            matchword = 'oxygenase'
        elif re.match('(.*)ogenase', word):
            matchword = 'ogenase'
        elif re.match('(.*)kinase', word):
            matchword = 'kinase'
        else:
            matchword = 'other'
    else:
        first_match = 'other'
        flag = 0
    if flag == 0:
        if first_match == 'other':
            flag_other = False
            for item in enzymes_dict['ase'][first_match]:
                for part in re.split('[ ,/()]',item):
                    if part.lower() == word: #.strip(',').strip('/').strip('(').strip(')'):
                        enzymes_list = enzymes_dict['ase'][first_match]
                        flag_other = True
                        break
                if flag_other == True:
                    break
        else:
            enzymes_list = enzymes_dict['ase'][first_match]
        
    else:
        if matchword == 'other':
            flag_other = False
            for item in enzymes_dict['ase'][first_match][matchword]:
                for part in re.split('[ ,/()]',item):
                    if part.lower() == word: #.strip(',').strip('/').strip('(').strip(')'):
                        enzymes_list = enzymes_dict['ase'][first_match][matchword]
                        flag_other = True
                        break
                if flag_other == True:
                    break
        else:
            enzymes_list = enzymes_dict['ase'][first_match][matchword]
    return enzymes_list
    
def part_match(text, keyword, search_list): # search using "re"--before and after
    entity = None
    words = text.strip().split() #re.split(' ', text.strip())
    index = 0
    text_after = ''
    text_before = ''
    for word in words:
        if word.find(keyword)>=0:
            # Search the words after keyword
            if word[-1] in end_word:
                entity = word[:-1].strip('s').strip('{').strip('[').strip('(').strip('(').strip('"')
                if entity[-1] in end_word:
                    entity = entity[:-1]
                text_after = ' '.join(words[index+1:]).strip()
            elif word[-4:] == 'ases':
                entity = word.strip('{').strip('[').strip('(').strip('"')
            else:
                entity = word.strip('{').strip('[').strip('(').strip('"')
                text_after = ' '.join(words[index+1:]).strip()
                pattern = re.compile(r'[(](.*?)[)]', re.I)
                after_bracket = pattern.match(text_after)
                if  after_bracket != None and after_bracket.group().strip('(').strip(')') in list_bracket:
                    entity = entity + ' ' + after_bracket.group()
                    text_after = text_after[after_bracket.end():].strip()
            #entity = word.strip(':').strip(',').strip('!').strip('?').strip('{').strip('[').strip('(').strip('"').strip('}').strip(']').strip(')').strip('"')

            # Search the word before keyword
            index_tmp = index
            for i in range(1, 5):
                flag = False
                if word[0] not in end_word and index - i >= 0 and words[index-i].strip('(').strip('[').strip('"') not in nonsense_list:
                    word_before = words[index-i].strip('(').strip('[').strip('"')
                    for item in search_list: 
                        parts = item.split()
                        for part in parts:
                            if part == word_before and not re.match('\d+', word_before):
                                entity = word_before + ' ' + entity
                                if entity[0] in end_word:
                                    entity = entity[1:]
                                del words[index-i]
                                index_tmp -= 1
                                flag = True
                                break
                        if flag == True:
                            break
                if flag == False:
                    break
            text_before = ' '.join(words[: index_tmp]).strip()
                
            break   
        index += 1
    text = text_before + ' ' + text_after
    
    global outside_list_num 
    outside_list_num += 1
    return entity, text, keyword

def search_ase(text, word):  
    search_list = enzyme_ase(word)
    flag = False # To mark if need to do part_match.
    entity = None
    if search_list:
        maxl = 0
        phrase = ''
        for item in search_list:
            pstart = text.find(item)
            #words = [i.strip('{').strip('[').strip('(').strip('(').strip('"') for i in text.split()]
            if pstart >= 0 and item.split()[0] in text.split(): # greedy 
                if len(item) > maxl:
                   entity = item
                   maxl = len(item)
                   phrase = text[pstart : len(text)]
        if entity != None:
            words = phrase.split()
            parts = entity.split()
            for i in range(len(parts)):
                if words[i] != parts[i]:
                    entity = None
                    return entity, text, 'EC numbers'
                    break
        if entity != None:
            global within_list_num
            start = text.find(entity)
            entity = (text[start:start + len(entity)])
            text = text[0 : start] + text[(start + len(entity) + 1) : len(text)] # There is least one space between two words.
            flag = True
            within_list_num += 1
    else: # we assume that enzymes are all the types (ending with existed 'ase') which can be find inside the kegg, if it not existed in the kegg, we do not use RE to check the text.
        entity = None
        flag = True
    if flag == False:
        return part_match(text, word, search_list)
    return entity, text, 'EC numbers'

def search_pattern(text, word, matchword):
    if matchword == 'other':
        entity = None
        if word.lower() in enzymes_dict[matchword]:
            entity = word
            start = text.lower().find(word)
            text = text[0:start].strip() + ' ' + text[start+len(entity)+1:len(sentence)].strip()
        return entity, text, 'EC numbers'
    if matchword != 'ase' and matchword != 'other':
        search_list = enzymes_dict[matchword]
        entity = None
        maxl = 0
        phrase = ''
        for item in search_list:
            pstart = text.find(item)
            if pstart>=0 and item.split()[0] in text.split(): # greedy
                if len(item) > maxl:
                   entity = item
                   maxl = len(item)
                   phrase = text[pstart : len(text)]
        if entity != None:
            words = phrase.split()
            parts = entity.split()
            for i in range(len(parts)):
                if words[i] != parts[i]:
                    entity = None
                    break
        if entity != None:
            global within_list_num
            start = text.find(entity)
            entity = text[start:start + len(entity)]
            text = text[0 : start] + text[(start + len(entity) + 1) : len(text)]
            within_list_num += 1
        return entity, text, 'EC numbers'
    else:
        return search_ase(text, word.rstrip('s'))
        
def annotation_dictionary():
    dict={'text':'', 
          'infons': {'identifier':'',
                     'type':'enzyme',
                     'annotator':'m.wang21@imperial.ac.uk',
                     'updated_at':time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())
                     },
          'id':'',
          'locations':[{'length':'',
                       'offset':''}]
                              
          }
    return dict
def find_all(entity, text):
    pos_list = []
    index = 0
    while index < len(text):
        1
    return pos_list
def annotate_corpus(id_num, ider, entity, offset, text, s_idx):
    if entity not in dict_entities:
        dict_entities[entity] = 0
    else:
        dict_entities[entity] += 1
    tmp = entity.replace('\\','\\\\\\').replace('(','\(').replace(')','\)').replace('[','\[').replace('{','\{').replace('+','\+')

    positions = [item.start() for item in re.finditer(tmp, text)]
    offset = offset + positions[dict_entities[entity]]
    #print (text[positions[dict_entities[entity]]:positions[dict_entities[entity]]+len(entity)])
    annotation_dict = annotation_dictionary()
    annotation_dict['text'] = entity
    annotation_dict['id'] = str(id_num)
    annotation_dict['infons']['identifier'] = ider
    annotation_dict['locations'][0]['length'] = len(entity)
    annotation_dict['locations'][0]['offset'] = offset

    return annotation_dict
'''
def abbreviation_replace(fulltext):
    
    return fulltext
'''

if __name__ == "__main__":
    start_time = time.time()

    enzymes_path = ("D:/BMR-DS/Project 1/Dataset/Classification.json")
    with open(enzymes_path, 'r', encoding = 'utf-8') as enzymes:
        enzymes_dict = json.load(enzymes)
    
    print_path = ("D:/BMR-DS/Project 1/Dataset/printnew.txt")
    w = open (print_path, 'a',encoding = 'utf-8')
    
    bracket_path = ("D:/BMR-DS/Project 1/Dataset/After_word.txt")
    with open(bracket_path, 'r') as bracket:
        list_bracket = bracket.readlines()
    bracket.close()
        
    # if there is a request to give the past time of a verb?
    notase_list = ['release','increase', 'database', 'base', 'disease','decrease', 'case']
    nonsense_list=['or', 'to', 'in', 'the', 'a', 'on', 'of', 'for', 'at', 'with','and','this','that','these','those','is','are']
    end_word = [')', '.', ';', '?', '!', '"', ':', ']', '}', '(', '[', '{',',','/']
    # Patterns used to rematch.
    pattern1 = re.compile(r'complex', re.I)
    pattern2 = re.compile(r'enzyme', re.I)
    pattern3 = re.compile(r'protein', re.I)
    pattern4 = re.compile(r'transporter', re.I)
    pattern5 = re.compile(r'Transferred')
    pattern6 = re.compile(r'doxin', re.I)
    pattern7 = re.compile(r'system', re.I)
    pattern8 = re.compile(r'cytochrome', re.I)
    pattern9 = re.compile(r'(.*)ase$|(.*)ases$')
    
    within_list_num = 0
    outside_list_num = 0
    file_num = 0
    for filename in os.listdir('D:/BMR-DS/Project 1/Dataset/Main_Text'):
        file_num += 1
        
        file_path = 'D:/BMR-DS/Project 1/Dataset/Main_Text/' + filename
        with open(file_path, 'r', encoding = 'utf-8') as file: # windows 10 has the 'gbk' codec problem without encoding = 'utf-8'
            fulltext = json.load(file)     
            
        # Abbreviation replacing    
        # fulltext = abbreviation_replace(fulltext)    
        
        #print(filename+":")
        w.write(filename)
        w.write('\n')
        
        within_list_num = 0
        outside_list_num = 0
        id_num = 0
        cont = 0 
        
        for text in fulltext['documents'][0]['passages']:
            paragraph = text['text'].split('.')
            annotation = text['annotations']
            offset = int(fulltext['documents'][0]['passages'][cont]['offset'])
            dict_entities = {}
            s_idx = 0 
            for sentence in paragraph:
                sentence = sentence.strip()
                words = re.split('[ ]', sentence)
                matchword = ''
                for word in words:
                    if pattern1.search(word):
                        matchword = 'complex'
                    elif pattern2.search(word):
                        matchword = 'enzyme'
                    elif pattern3.search(word):
                        matchword = 'protein'
                    elif pattern4.search(word):
                        matchword = 'transporter'
                    elif pattern5.search(word):
                        matchword = 'Transferred'
                    elif pattern6.search(word):
                        matchword = 'doxin'
                    elif pattern7.search(word):
                        matchword = 'system'
                    elif pattern8.search(word):
                        matchword = 'cytochrome'
                    elif pattern9.search(word) and pattern9.search(word).group().rstrip('s').lower() not in notase_list:
                        matchword = 'ase'
                    else:
                        matchword = 'other'
                    (entity, content, ider) = search_pattern(sentence, word, matchword)
                    if entity != None:
                        w.write(entity + ' ' + ider)
                        w.write('\n')
                        id_num += 1
                        annotation.append(annotate_corpus(id_num, ider, entity, offset, text['text'], s_idx))
                        sentence = content
                s_idx += 1         
            fulltext['documents'][0]['passages'][cont]['annotations'] = annotation
            cont += 1
        '''
        if id_num == 0:
           print (' There is no enzyme found in '+filename+'.')
        else:
            print (' Totally, '+str(id_num)+' enzymes has been found in '+filename+'.')
            s = ' (in the list, out of list) = ('+str(within_list_num)+', '+str(outside_list_num)+')'
            print(s)
        '''
        file.close()
        writein_path = "D:/BMR-DS/Project 1/Processing/" + filename.rstrip('.json') + '_' + 'annotation.json'
        writein = open(writein_path, 'w')
        json.dump(fulltext, writein, indent = 4)
        writein.close()
    
    enzymes.close()
    writein.close()
    file.close()
    w.close()

    
    print('The total time of running the folder is: ' + str(time.time() - start_time))
    # print('The average time of each file is: ' + str((timef_end-timef_start)/ file_num))
    


'''

#import numpy as np
#import pandas as pd
import json
import re
import os
import time

def enzyme_ase(word):
    matchword = ''
    first_match = ''
    enzymes_list = []
    flag = 1 # If flag==0, search_list only has first_match.
    if re.match('(.*)synthase', word):
        first_match = 'synthase'
        flag = 0
    elif re.match('(.*)lyase', word):
        first_match = 'lyase'
        flag = 0
    elif re.match('(.*)ligase', word):
        first_match = 'ligase'
        flag = 0
    elif re.match('(.*)ease', word):
        first_match = 'ease'
        if re.match('(.*)permease', word):
            matchword = 'permease'
        elif re.match('(.*)protease', word):
            matchword = 'protease'
        else:
            matchword = 'other'
    elif re.match('(.*)lase', word):
        first_match = 'lase'
        if re.match('(.*)methylase', word):
            matchword = 'methylase'
        elif re.match('(.*)horylase', word):
            matchword = 'horylase'
        elif re.match('(.*)cyclase', word):
            matchword = 'cyclase'
        elif re.match('(.*)hydrolase', word):
            matchword = 'hydrolase'
        elif re.match('(.*)oxylase', word):
            matchword = 'oxylase'
        else:
            matchword = 'other'
    elif re.match('(.*)rase', word):
        first_match = 'rase'
        if re.match('(.*)transferase', word):
            matchword = 'transferase'
        elif re.match('(.*)saturase', word):
            matchword = 'saturase'
        else:
            matchword = 'other'
    elif re.match('(.*)tase', word):
        first_match = 'tase'
        if re.match('(.*)synthetase', word):
            matchword = 'synthetase'
        elif re.match('(.*)reductase(.*)', word):
            matchword = 'reductase'
        else:
            matchword = 'other'
    elif re.match('(.*)idase', word):
        first_match = 'idase'
        if re.match('(.*)oxidase', word):
            matchword = 'oxidase'
        elif re.match('(.*)peptidase', word):
            matchword = 'peptidase'
        else:
            matchword = 'other'
    elif re.match('(.*)nase', word):
        first_match = 'nase'
        if re.match('(.*)proteinase', word):
            matchword = 'proteinase'
        elif re.match('(.*)oxygenase', word):
            matchword = 'oxygenase'
        elif re.match('(.*)ogenase', word):
            matchword = 'ogenase'
        elif re.match('(.*)kinase', word):
            matchword = 'kinase'
        else:
            matchword = 'other'
    else:
        first_match = 'other'
        flag = 0
    if flag == 0:
        if first_match == 'other':
            flag_other = False
            for item in enzymes_dict['ase'][first_match]:
                for part in re.split('[ ,/()]',item):
                    if part.lower() == word: #.strip(',').strip('/').strip('(').strip(')'):
                        enzymes_list = enzymes_dict['ase'][first_match]
                        flag_other = True
                        break
                if flag_other == True:
                    break
        else:
            enzymes_list = enzymes_dict['ase'][first_match]
        
    else:
        if matchword == 'other':
            flag_other = False
            for item in enzymes_dict['ase'][first_match][matchword]:
                for part in re.split('[ ,/()]',item):
                    if part.lower() == word: #.strip(',').strip('/').strip('(').strip(')'):
                        enzymes_list = enzymes_dict['ase'][first_match][matchword]
                        flag_other = True
                        break
                if flag_other == True:
                    break
        else:
            enzymes_list = enzymes_dict['ase'][first_match][matchword]
    return enzymes_list
    
def part_match(text, keyword, search_list): # search using "re"--before and after
    entity = None
    words = text.strip().split() #re.split(' ', text.strip())
    index = 0
    text_after = ''
    text_before = ''
    for word in words:
        if word.find(keyword)>=0:
            # Search the words after keyword
            if word[-1] in end_word:
                entity = word[:-1].strip('s').strip('{').strip('[').strip('(').strip('(').strip('"')
                if entity[-1] in end_word:
                    entity = entity[:-1]
                text_after = ' '.join(words[index+1:]).strip()
            elif word[-4:] == 'ases':
                entity = word.strip('{').strip('[').strip('(').strip('"')
            else:
                entity = word.strip('{').strip('[').strip('(').strip('"')
                text_after = ' '.join(words[index+1:]).strip()
                pattern = re.compile(r'[(](.*?)[)]', re.I)
                after_bracket = pattern.match(text_after)
                if  after_bracket != None and after_bracket.group().strip('(').strip(')') in list_bracket:
                    entity = entity + ' ' + after_bracket.group()
                    text_after = text_after[after_bracket.end():].strip()
            #entity = word.strip(':').strip(',').strip('!').strip('?').strip('{').strip('[').strip('(').strip('"').strip('}').strip(']').strip(')').strip('"')

            # Search the word before keyword
            index_tmp = index
            for i in range(1, 5):
                flag = False
                if word[0] not in end_word and index - i >= 0 and words[index-i].strip('(').strip('[').strip('"') not in nonsense_list:
                    word_before = words[index-i].strip('(').strip('[').strip('"')
                    for item in search_list: 
                        parts = item.split()
                        for part in parts:
                            if part == word_before and not re.match('\d+', word_before):
                                entity = word_before + ' ' + entity
                                if entity[0] in end_word:
                                    entity = entity[1:]
                                del words[index-i]
                                index_tmp -= 1
                                flag = True
                                break
                        if flag == True:
                            break
                if flag == False:
                    break
            text_before = ' '.join(words[: index_tmp]).strip()
                
            break   
        index += 1
    text = text_after #text_before + ' ' + text_after
    
    global outside_list_num 
    outside_list_num += 1
    return entity, text, keyword

def search_ase(text, word):  
    search_list = enzyme_ase(word)
    flag = False # To mark if need to do part_match.
    entity = None
    if search_list:
        maxl = 0
        phrase = ''
        for item in search_list:
            pstart = text.lower().find(item.lower())
            #words = [i.strip('{').strip('[').strip('(').strip('(').strip('"') for i in text.split()]
            if pstart >= 0 and item.split()[0].lower() in text.lower().split(): # greedy 
                if len(item) > maxl:
                   maxl = len(item)
                   entity = text[pstart:pstart+maxl]
                   phrase = text[pstart : len(text)]
        if entity != None:
            words = phrase.split()
            parts = entity.split()
            for i in range(len(parts)):
                if words[i] != parts[i]:
                    entity = None
                    return entity, text, 'EC numbers'
                    break
        if entity != None:
            global within_list_num
            start = text.lower().find(entity.lower())
            #entity = text[start:start + len(entity)]
            text = text[(start + len(entity) + 1) :].strip()#text[0 : start].strip() + ' ' + text[(start + len(entity) + 1) : len(text)].strip() # There is least one space between two words.
            flag = True
            within_list_num += 1
    else: # we assume that enzymes are all the types (ending with existed 'ase') which can be find inside the kegg, if it not existed in the kegg, we do not use RE to check the text.
        entity = None
        flag = True
    if flag == False:
        return part_match(text, word, search_list)
    return entity, text, 'EC numbers'

def search_pattern(text, word, matchword):
    if matchword == 'other':
        entity = None
        if word.lower() in enzymes_dict[matchword]:
            entity = word
            start = text.lower().find(word.lower())
            text = text[start+len(entity)+1:len(sentence)].strip()#text[0:start].strip() + ' ' + text[start+len(entity)+1:len(sentence)].strip()
        return entity, text, 'EC numbers'
    if matchword != 'ase' and matchword != 'other':
        search_list = enzymes_dict[matchword]
        entity = None
        maxl = 0
        phrase = ''
        for item in search_list:
            pstart = text.lower().find(item.lower())
            if pstart>=0 and item.split()[0].lower() in text.lower().split(): # greedy
                if len(item) > maxl:
                   maxl = len(item)
                   entity = text[pstart:pstart+maxl]
                   phrase = text[pstart : len(text)]
        if entity != None:
            words = phrase.split()
            parts = entity.split()
            for i in range(len(parts)):
                if words[i] != parts[i]:
                    entity = None
                    break
        if entity != None:
            global within_list_num
            within_list_num += 1
            start = text.lower().find(entity.lower())
            #entity = text[start:start + len(entity)]
            text = text[(start + len(entity) + 1) : len(text)].strip()#text[0 : start].strip()+ ' ' + text[(start + len(entity) + 1) : len(text)].strip()
            
        return entity, text, 'EC numbers'
    else:
        return search_ase(text, word.rstrip('s'))
        
def annotation_dictionary():
    dict={'text':'', 
                 'infons': {'identifier':'',
                            'type':'enzyme',
                            'annotator':'m.wang21@imperial.ac.uk',
                            'updated_at':time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())
                            },
                 'id':'',
                 'locations':[{'length':'',
                              'offset':''
                              }]
                              
                }
    return dict
def find_all(entity, text):
    pos_list = []
    index = 0
    while index < len(text):
        1
    return pos_list
def annotate_corpus(id_num, ider, entity, offset, text, s_idx):
    if entity not in dict_entities:
        dict_entities[entity] = 0
    else:
        dict_entities[entity] += 1
    tmp = entity.replace('\\','\\\\\\').replace('(','\(').replace(')','\)').replace('[','\[').replace('{','\{').replace('+','\+')
    print (entity)
    print (tmp)
    positions = [item.start() for item in re.finditer(tmp, text.replace('\n',' '))]
    print (positions)
    offset = offset + positions[dict_entities[entity]]
    #print (text[positions[dict_entities[entity]]:positions[dict_entities[entity]]+len(entity)])
    annotation_dict = annotation_dictionary()
    annotation_dict['text'] = entity
    annotation_dict['id'] = str(id_num)
    annotation_dict['infons']['identifier'] = ider
    annotation_dict['locations'][0]['length'] = len(entity)
    annotation_dict['locations'][0]['offset'] = offset

    return annotation_dict
'''
def abbreviation_replace(fulltext):
    
    return fulltext
'''

if __name__ == "__main__":
    start_time = time.time()

    enzymes_path = ("Auto-CORPus/Classification.json")
    with open(enzymes_path, 'r', encoding = 'utf-8') as enzymes:
        enzymes_dict = json.load(enzymes)
    
    print_path = ("Auto-CORPus/printnew-MWAS.txt")
    w = open (print_path, 'a',encoding = 'utf-8')
    
    bracket_path = ("Auto-CORPus/After_word.txt")
    with open(bracket_path, 'r') as bracket:
        list_bracket = bracket.readlines()
    bracket.close()
        
    # if there is a request to give the past time of a verb?
    notase_list = ['release','increase', 'database', 'base', 'disease','decrease', 'case']
    nonsense_list=['or', 'to', 'in', 'the', 'a', 'on', 'of', 'for', 'at', 'with','and','this','that','these','those','is','are']
    end_word = [')', '.', ';', '?', '!', '"', ':', ']', '}', '(', '[', '{',',','/']
    # Patterns used to rematch.
    pattern1 = re.compile(r'complex', re.I)
    pattern2 = re.compile(r'enzyme', re.I)
    pattern3 = re.compile(r'protein', re.I)
    pattern4 = re.compile(r'transporter', re.I)
    pattern5 = re.compile(r'Transferred')
    pattern6 = re.compile(r'doxin', re.I)
    pattern7 = re.compile(r'system', re.I)
    pattern8 = re.compile(r'cytochrome', re.I)
    pattern9 = re.compile(r'(.*)ase$|(.*)ases$')
    
    within_list_num = 0
    outside_list_num = 0
    file_num = 0
    enzyme_num = 0
    
    for filename in os.listdir('Auto-CORPus/output-MWAS/Main_Text'):
        
        
        if re.match('(.*).json$', filename):
            file_path = 'Auto-CORPus/output-MWAS/Main_Text/' + filename
            with open(file_path, 'r', encoding = 'utf-8') as file: # windows 10 has the 'gbk' codec problem without encoding = 'utf-8'
                fulltext = json.load(file)    
            file_num += 1
        else:
            continue
        # Abbreviation replacing    
        # fulltext = abbreviation_replace(fulltext)    
        
        print(filename+":")
        w.write(filename)
        w.write('\n')
        
        within_list_num = 0
        outside_list_num = 0
        id_num = 0
        cont = 0 
        
        for text in fulltext['documents'][0]['passages']:
            paragraph = text['text'].split('.')
            annotation = text['annotations']
            offset = int(fulltext['documents'][0]['passages'][cont]['offset'])
            dict_entities = {}
            s_idx = 0 
            for sentence in paragraph:
                sentence = sentence.strip()
                words = re.split('[ ]', sentence)
                matchword = ''
                for word in words:
                    if pattern1.search(word):
                        matchword = 'complex'
                    elif pattern2.search(word):
                        matchword = 'enzyme'
                    elif pattern3.search(word):
                        matchword = 'protein'
                    elif pattern4.search(word):
                        matchword = 'transporter'
                    elif pattern5.search(word):
                        matchword = 'Transferred'
                    elif pattern6.search(word):
                        matchword = 'doxin'
                    elif pattern7.search(word):
                        matchword = 'system'
                    elif pattern8.search(word):
                        matchword = 'cytochrome'
                    elif pattern9.search(word) and pattern9.search(word).group().rstrip('s').lower() not in notase_list:
                        matchword = 'ase'
                    else:
                        matchword = 'other'
                    (entity, content, ider) = search_pattern(sentence, word, matchword)
                    if entity != None:
                        w.write(entity + ' ' + ider)
                        w.write('\n')
                        id_num += 1
                        annotation.append(annotate_corpus(id_num, ider, entity, offset, text['text'], s_idx))
                        sentence = content
                s_idx += 1         
            fulltext['documents'][0]['passages'][cont]['annotations'] = annotation
            cont += 1
        enzyme_num += id_num
        '''
        if id_num == 0:
           print (' There is no enzyme found in '+filename+'.')
        else:
            print (' Totally, '+str(id_num)+' enzymes has been found in '+filename+'.')
            s = ' (in the list, out of list) = ('+str(within_list_num)+', '+str(outside_list_num)+')'
            print(s)
        '''
        file.close()
        writein_path = "Auto-CORPus/Annotated_Corpus/MWAS/" + filename.rstrip('.json') + '_' + 'annotated.json'
        writein = open(writein_path, 'w')
        json.dump(fulltext, writein, indent = 4)
        writein.close()
    
    enzymes.close()
    w.close()

    
    print('The total time of running {} files is: '.format(file_num) + str(time.time() - start_time) + '.')
    print('{} enzymes have been found in this folder.'.format(enzyme_num))
    # print('The average time of each file is: ' + str((timef_end-timef_start)/ file_num))
    
'''