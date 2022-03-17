#### -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 22:00:51 2022

@author: one
"""

import os 
import json
import re
import time

def enzyme_ase_list(word):
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
        elif re.match('(.*)reductase', word):
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

def search_ase(text, word):  
    global isEC
    global oritext
    global pri
    search_list = enzyme_ase_list(word.lower())
    flag = False # To mark if need to do part_match.
    entity = None
    if search_list:
        maxl = 0
        phrase = ''
        words = []
        for i in text.lower().split():
            if i == '/':
                continue
            words.append(i[i.find('/')+1:])
        for item in search_list:
            pstart = text.lower().find(item.lower())
            if pstart >= 0 and (pstart == 0 or text[pstart-1] == ' ' or text[pstart-1] == ':'): #text.lower().split(): # greedy 
                if len(item) > maxl:
                    maxl = len(item)
                    entity = text[pstart : pstart+maxl]
                    phrase = text[pstart :]
        if entity != None:
            words = phrase.split()
            parts = entity.split()
            for i in range(len(parts)):
                if i == len(parts) -1 and words[i][-1] in [',',';','!','"', "'", '?']:
                    words[i] = words[i][:-1]
                elif i == len(parts) -1 and words[i].find('/') >= 0:
                    words[i] = words[i][:words[i].find('/')]
                elif i == len(parts) -1 and words[i][-4:] == 'ases':
                    entity = entity + 's'
                    parts[i] = words[i]
                if words[i] != parts[i]:
                    entity = None
                    return entity, text, 'EC numbers'
                    break
        if entity != None:
            start = text.lower().find(entity.lower())
            end = start + len(entity)
            words = text.lower().split()
            if start - 1 >= 0 and oritext[start-1] == '-':
                index = 0
                for index in range(len(words)):
                    lEntity = len(entity.split())
                    #if re.search(entity.lower(), ' '.join(words[index:(index+lEntity)])):
                    if entity.lower() == ' '.join(words[index:(index+lEntity)]).strip(';').strip(':').strip(','):
                        break
                #index = words.index(entity.split()[-1].lower())
                #index = index - len(entity.split()) + 1
                if words[index-1][-1] != ')' or ((words[index-1][-1] == ')' and re.search('[(].*[)]',words[index-1]))):# or re.match('[(](.*)[)]$|[(](.*)[)],$', ' '.join(words[index-1:(index+lEntity)])):
                    start -= (len(words[index-1])+1)
                    entity = text[start : end]
                    if entity[0] == '(' and not re.match('[(].*[)]',entity):#')' not in entity:
                        entity = entity[1:]
                        start += 1
                    end = start+len(entity)#+len(words[index-1])+1
                    isEC = False
            if re.search('RNA polymerase III', entity):
                print(entity)
            if end < len(text) and oritext[end] == '-':
                for index in range(len(words)):
                    lEntity = len(entity.split())
                    if entity.lower() == ' '.join(words[index:(index+lEntity)]).strip(';').strip(':').strip(',') :
                        break
                index = index + lEntity -1
                #index = words.index(entity.split()[-1].lower())
                if index+1 < len(words) and (re.match('\d|(.*)DNA$|(.*)RNA$', words[index+1]) or words[index+1] in pri):
                    end = start+len(entity)+len(words[index+1])+1
                    entity = text[start : end]
                    entity = entity.strip(',').strip(']').strip(')').strip(',').strip(';')
                    isEC = False
            text = text[(start + len(entity) + 1):].strip()
            oritext = oritext[(start + len(entity) + 1):].strip()
    else:
        entity = None
    if entity == 'RNA polymerase III':
        print(entity)
    return entity, text, 'EC numbers'


def search_pattern(text, word, matchword):
    global oritext
    global isEC 
    isEC = True
    entity = None
    if word.find('/') >= 0:
        word = word[word.find('/')+1:]
    if matchword == 'other':
        if word.lower() in enzymes_dict[matchword]:
            entity = word
            start = text.lower().find(word.lower())
            text = text[start+len(entity)+1:].strip()#text[0:start].strip() + ' ' + text[start+len(entity)+1:len(sentence)].strip()
        return entity, text, 'EC numbers'
    if matchword != 'ase' and matchword != 'other':
        search_list = enzymes_dict[matchword]
        maxl = 0
        phrase = ''
        words = []
        for i in text.lower().split():
            if i == '/':
                continue
            words.append(i[i.find('/')+1:])
        for item in search_list:
            pstart = text.lower().find(item.lower())
            if  pstart >= 0 and (pstart == 0 or text[pstart-1] == ' ' or text[pstart-1] == ':'): # greedy
                if len(item) > maxl:
                    maxl = len(item)
                    entity = text[pstart : pstart+maxl]
                    phrase = text[pstart :]
        if entity != None:
            words = phrase.split()
            parts = entity.split()
            for i in range(len(parts)):
                if i == len(parts) -1 and words[i][-1] in [',',';','!','"', "'", '?']:
                    words[i] = words[i][:-1]
                elif i == len(parts) -1 and words[i].find('/') >= 0:
                    words[i] = words[i][:words[i].find('/')]
                if words[i] != parts[i]:
                    entity = None
                    break
        if entity != None:
            start = text.lower().find(entity.lower())
            end = start + len(entity)
            words = text.lower().split()
            if start - 1 >= 0 and oritext[start-1] == '-':
                index = 0
                for index in range(len(words)):
                    lEntity = len(entity.split())
                    if re.search(entity.lower(), ' '.join(words[index:(index+lEntity)])):
                        break
                #index = words.index(entity.split()[-1].lower())
                #index = index - len(entity.split()) + 1
                if words[index-1][-1] != ')' or ((words[index-1][-1] == ')' and re.search('[(].*[)]',words[index-1]))):
                    start -= (len(words[index-1])+1)
                    entity = text[start : end]
                    if entity[0] == '(' and not re.match('[(].*[)]',entity):#')' not in entity:#in end_word:
                        entity = entity[1:]
                        start += 1
                    end = start+len(entity)+len(words[index-1])+1
                    isEC = False
            if end < len(text) and oritext[end] == '-':
                for index in range(len(words)):
                    lEntity = len(entity.split())
                    if re.search(entity.lower(), ' '.join(words[index:(index+lEntity)])):
                        break
                index = index + lEntity -1
                #index = words.index(entity.split()[-1].lower())
                if index+1 < len(words) and (re.match('\d|(.*)DNA$|(.*)RNA$', words[index+1]) or words[index+1] in pri):
                    end = start+len(entity)+len(words[index+1])+1
                    entity = text[start : end]
                    entity = entity.strip(',').strip(']').strip(')').strip(',').strip(';')
                    isEC = False
            text = text[(start + len(entity) + 1):].strip()
            oritext = oritext[(start + len(entity) + 1):].strip()
        return entity, text, 'EC numbers'
    else:
        return search_ase(text, word.rstrip('s'))


def search_enzyme(word):
    matchword = ''
    if pattern1.search(word):
        matchword = 'complex'
    elif pattern9.search(word) and pattern9.search(word).group().rstrip('s').lower() not in notase_list:
        matchword = 'ase'    
    elif pattern6.search(word):
        matchword = 'doxin'
    elif pattern2.search(word):
        matchword = 'enzyme'
    elif pattern4.search(word):
        matchword = 'transporter'
    elif pattern5.search(word):
        matchword = 'Transferred'
    elif pattern7.search(word):
        matchword = 'system'
    elif pattern8.search(word):
        matchword = 'cytochrome'    
    elif pattern3.search(word):
        matchword = 'protein'
    else:
        matchword = 'other'
    return matchword

def annotation_dictionary():
    dict={'text':'', 
          'infons': {'identifier':'',
                     'type':'enzyme',
                     'annotator':'m.wang21@fuzzysearch',
                     'updated_at':time.strftime("%Y-%m-%dT%H:%M:%SZ",time.localtime())
                    },
          'id':'',
          'locations':[{'length':'',
                        'offset':''
                      }]
                              
          }
    return dict

def annotate_corpus(ider, entity, originaloffset, text, originaltext, section, annotation):
    global ec_num
    global kw_num
    global id_num
    global isEC
    global num
    '''
    if entity not in dict_entities:
        dict_entities[entity] = 0
    else:
        dict_entities[entity] += 1
    tmp = entity.replace('\\','\\\\\\').replace('(','\(').replace(')','\)').replace('[','\[').replace('{','\{').replace('+','\+').replace('*', '\*')
    positions = [item.start() for item in re.finditer(tmp, text.replace('\n',' '))]
    position = positions[dict_entities[entity]]
    offset = originaloffset + position
    '''
    
    position = text.lower().find(entity.lower())
    offset = originaloffset + position
    if entity == 'serine/threonine protein kinase':
        print(text)
        print(entity)
        print(position)
        print(offset)
    length = len(entity)
    print(originaltext[position : (position + length)])
    position, offset, length = double_check(originaltext, position, position + length, originaloffset)
    #print(originaltext[position : (position + length)])
    flag = False
    if annotation == []:
        flag == True
    for item in annotation:
        ioffset = item['locations'][0]['offset']
        ilength =  item['locations'][0]['length']
        if offset == ioffset:
            if length > ilength:
                annotation.remove(item)
                id_num -= 1
                flag = True
            elif ilength == length:
                flag = False
                break
        elif offset > ioffset and originaltext[position - 1] == '-':
            if offset + length < ioffset + ilength:
                flag = False
                break
            elif offset + length == ioffset + ilength:
                length = length + offset - item['locations'][0]['offset']
                offset = item['locations'][0]['offset']
                position = item['locations'][0]['offset'] - originaloffset
                annotation.remove(item)
                id_num -= 1
                isEC = False
                flag = True
                break
            elif offset < ioffset + ilength and offset + length > ioffset + ilength:
                length = length + offset - item['locations'][0]['offset']
                offset = item['locations'][0]['offset']
                position = item['locations'][0]['offset'] - originaloffset
                annotation.remove(item)
                id_num -= 1
                isEC = False
                flag = True
                break
        elif offset > ioffset and offset + length == ioffset + ilength:
            flag = False
            break
        elif offset < ioffset and offset + length == ioffset + ilength:
            if text[ioffset-originaloffset-1] != '/':
                annotation.remove(item)
                id_num -= 1
                flag = True
                break
            else:
                flag = False
                break
        else:
            flag = True
    if flag == True:
        id_num += 1
        num += 1
        annotation_dict = annotation_dictionary()
        annotation_dict['id'] = str(id_num)
        if isEC == True:
            annotation_dict['infons']['identifier'] = ec_number_list[entity.rstrip('s').lower()]
        else:
            annotation_dict['infons']['identifier'] = 'fuzzy_search'
        entity = originaltext[position : (position+length)]
        annotation_dict['text'] = entity
        annotation_dict['locations'][0]['length'] = len(entity)
        annotation_dict['locations'][0]['offset'] = offset
        #print(entity)
        return annotation_dict, annotation
    else:
        return None, annotation

def double_check(originaltext, oriposition, end, originaloffset):
    s = originaltext[oriposition]
    if s in ['"', "'", '/','‘', '“', '[', '{']:
        return oriposition+1, originaloffset, end-oriposition
    elif s == '(' and not re.match(r'[(](.*?)[)]', originaltext[oriposition : end]):
        return oriposition+1, originaloffset, end-oriposition
    #if originaltext[oriposition] in ['"', "'", '(', '[', '{', '/','‘', '“']:
    #    return oriposition, originaloffset, end-oriposition
    tmp = oriposition - 1
    position = oriposition
    flag = True
    global isEC
    while tmp > 0:
        if originaltext[tmp] == ' ':
            position = tmp + 1
            break
        elif originaltext[tmp] in ['"', "'", '/','(', '[', '{', '‘', '“']:
            position = tmp + 1
            flag = False
            return position, originaloffset, end-oriposition
        tmp -= 1
    # check the word before entity we found
    if tmp > 0 and len(originaltext.strip()) > 0:
        #print(originaltext)
        befword = originaltext[:tmp].strip().split()[-1]    
        f = False
        if befword[-1] not in ['.', ';', ':', ')', ']', '}', '!', '?', ',', '’', '”'] and flag:
            parts = befword.split('-')
            for part in parts:
                if part not in nonsense_list and not re.match('\d+', part):
                    for item in ec_number_list:
                        if part in item.split():
                            f = True
                            break
                    if f == True:
                        break
        if f == True:
            befword = befword.lstrip('"').lstrip("'")
            if befword[0] == '(' and not re.match(r'[(](.*?)[)]', befword):
                befword = befword.lstrip('(')
            elif befword[0] == '[' and not re.match(r'[[](.*?)[]]', befword):
                befword = befword.lstrip('[')
            position = position - len(befword) - 1
    if position != oriposition:
        isEC = False
    offset = position + originaloffset
    length = end - position
    return position, offset, length

if __name__ == "__main__":
    
    start_time = time.time()
    
    enzymes_path = ("Auto-CORPus/Classification_new.json")
    with open(enzymes_path, 'r', encoding = 'utf-8') as enzymes:
        enzymes_dict = json.load(enzymes)
        enzymes.close()
    
    ec_number_path = ("Auto-CORPus/KEGG_EC_new.json")
    with open(ec_number_path, 'r', encoding = 'utf-8') as ec_number:
        ec_number_list = json.load(ec_number)
    ec_number.close()
    
    notase_list = ['release','increase', 'database', 'base', 'disease','decrease', 'case']
    nonsense_list=['as', 'between', 'specific', 'new', 'old', 'an', 'or', 'to', 'in', 'the', 'a', 'on', 'of', 'for', 'at', 'with','and','this','that','these','those','is','are', 'enzyme', 'enzymes']
    end_word = [')', '.', ';', '?', '!', '"', ':', ']', '}', '(', '[', '{',',','/',]
    pri = ['sigma','alpha', 'beta', 'delta', 'gamma', 'zeta', 'epsilon', 'eta', 'theta', 'lota', 'kappa', 'lambda', 'mu','nu','xi','omicron','pi','rho','tau','upsilon','phi','chi','psi','omega']
    for i in range(945, 970):
        pri.append(chr(i))
    for i in range(913, 938):
        pri.append(chr(i))
    pattern1 = re.compile(r'complex', re.I)
    pattern2 = re.compile(r'enzyme', re.I)
    pattern3 = re.compile(r'protein', re.I)
    pattern4 = re.compile(r'transporter', re.I)
    pattern5 = re.compile(r'Transferred')
    pattern6 = re.compile(r'doxin', re.I)
    pattern7 = re.compile(r'system', re.I)
    pattern8 = re.compile(r'cytochrome', re.I)
    pattern9 = re.compile(r'(.*)ase$|(.*)ases$', re.I)
    
    enzyme_num = 0
    num = 0
     
    for filename in os.listdir('Auto-CORPus/Annotated_Corpus/proteomics-Abbre'):
        if re.match('(.*).json$', filename):
            file_path = 'Auto-CORPus/Annotated_Corpus/proteomics-Abbre/' + filename
            print(filename)
            with open(file_path, 'r', encoding = 'utf-8') as file: # windows 10 has the 'gbk' codec problem without encoding = 'utf-8'
                fulltext = json.load(file)    
            file.close() 
        else:
            continue
        flag = False
        cont = 0    
        id = 0
        for text in fulltext['documents'][0]['passages']:
            hyphen_Location = re.finditer('-', text['text'].replace('\n', ' '))
            paragraph = text['text'].split('.')

            annotation = text['annotations']
            id_num = len(annotation)
            offset = int(fulltext['documents'][0]['passages'][cont]['offset'])
            section = fulltext['documents'][0]['passages'][cont]['infons']['iao_name_1']
            dict_entities = {}
            dict_entities_abbre = {}
            curoffset = offset
            nextstart = 0
            for sentence in paragraph:
                oritext = sentence.strip()
                sentence = sentence.strip().replace('-', ' ')
                words = re.split('[ ]', sentence)
                for word in words:
                    matchword = search_enzyme(word.lower())
                    if matchword != '':
                        (entity, content, ider) = search_pattern(sentence, word, matchword)
                        if entity != None and entity != '':
                            item, annotation = annotate_corpus(ider, entity, curoffset, text['text'][nextstart:].replace('-', ' '), text['text'][nextstart:], section, annotation)
                            curoffset = curoffset + text['text'][nextstart:].replace('-', ' ').find(entity) + len(entity)
                            nextstart = nextstart + text['text'][nextstart:].replace('-', ' ').find(entity) + len(entity)
                            if item != None:
                                annotation.append(item)
                            sentence = content
                            
            for i in range(0, len(annotation)):
                id += 1
                annotation[i]['id'] = id
            
            fulltext['documents'][0]['passages'][cont]['annotations'] = annotation
            cont += 1
        
        enzyme_num += id
        if id > 0:
            writein_path = "Auto-CORPus/Annotated_Corpus/proteomics-greek/" + filename #.replace('_bioc_annotated.json', '_hyphen.json')
            writein = open(writein_path, 'w')
            json.dump(fulltext, writein, indent = 4)
            writein.close()      
    print('Time is {}.'.format(time.time()-start_time))
    print('{} enzymes have been found in the proteomics corpus.'.format(enzyme_num))
    print('{} enzymes have been found through fuzzy search.'.format(num))