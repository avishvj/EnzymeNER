# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:52:40 2022

This script is used to segment the text to sentence which is suitable for input of NER models

@author: Meiqi
"""

import spacy
import json
import pandas as pd
import time
import re
import os
import operator

def corpus_output(file_id, sec_id, subsec_id, sent_id, text):
    output_dict = {'corpus': [], 'section': [], 'text': []}
    
    subsec_str = str(subsec_id).zfill(2)+str(sent_id).zfill(3)
    output_dict['corpus'].append(file_id)
    output_dict['section'].append(sec_id + subsec_str)
    output_dict['text'].append(text)
    
    return pd.DataFrame(output_dict)

def annotation_output(file_id, sec_id, subsec_id, sent_id, start, end, entity):
    output_dict = {'corpus': [], 'section': [], 'start': [], 'end': [], 'enzyme':[]}
    
    subsec_str = str(subsec_id).zfill(2)+str(sent_id).zfill(3)
    output_dict['corpus'].append(file_id)
    output_dict['section'].append(sec_id + subsec_str)
    output_dict['start'].append(start)
    output_dict['end'].append(end)
    output_dict['enzyme'].append(entity)
    
    return pd.DataFrame(output_dict)
def files(fileName):
    fdlist = ['GWAS', 'MWAS', 'Microbiome', 'proteomics']
    for fdname in fdlist:
        folderPath = "D:/BMR-DS/Project 1/Dataset/{}".format(fdname)
        flist = os.listdir(folderPath)
        if fileName in flist:
            abbrevPath = 'D:/BMR-DS/Project 1/Dataset/{}-Abbreviations'.format(fdname)
            print(fdname)
            return folderPath, abbrevPath

if __name__ == '__main__':
    
    # segment the sentences
    nlp = spacy.load("en_core_web_sm")
    doc = nlp("This is a sentence. This is another sentence.")
    assert doc.has_annotation("SENT_START")
    # output files
    outputText_path = "D:/BMR-DS/Project 1/Dataset/Annalysis.txt"
    outputAnnot_path = "D:/BMR-DS/Project 1/Dataset/AnnotAnnalysis.txt"
    df_text = pd.DataFrame()
    df_annot = pd.DataFrame()
    
    sectionType = {'methods section': 'M',
                   'results section': 'R',
                   'discussion section': 'D',
                   'textual abstract section': 'A',
                   'introduction section': 'I',
                   'references section': 'F' 
                   }
    file_num = 0
    startTime = time.time()
    annonum = 0
    '''
    fmpath = 'D:/BMR-DS/Project 1/Dataset/ManualAnnotationFiles.txt'
    fm = open(fmpath, 'r')
    fmnames = []
    for i in fm.readlines():
        fmnames.append(i.replace('\n', '.json'))
    fm.close()
    print(fmnames)
    '''
    enzdict = {}
    fdlist = ['GWAS', 'MWAS', 'Microbiome', 'proteomics']
    for fdname in fdlist:
        folderPath = "D:/BMR-DS/Project 1/Dataset/{}".format(fdname)
        fmnames = os.listdir(folderPath)
        
        abbrevPath = 'D:/BMR-DS/Project 1/Dataset/{}-Abbreviations'.format(fdname)
        
        for fileName in fmnames:
            if re.match('(.*).json$', fileName):
                fileName = fileName.replace('_bioc_annotated.json', '_bioc_annotated.json')
            else:
                continue
            if re.match('(.*).json$', fileName):
                filePath = '/'.join([folderPath, fileName])
                with open(filePath, 'r', encoding = 'utf-8') as file:
                    mainText = json.load(file)
                file.close()
                abbrePath = '/'.join([abbrevPath, fileName.replace('bioc_annotated', 'abbreviations')])
                with open(abbrePath, 'r', encoding = 'utf-8') as abbrefile:
                    abbreText = json.load(abbrefile)
                abbrefile.close()
        
            else:
                continue
            
            file_num += 1
            fileID = fileName.split('_')[0]
            print('\r'+'Processing {} ({}/{})'.format(fileID, str(file_num), str(len(fmnames))), flush = True)
            
            sectionDict = {'M': -1, 'R': -1, 'D': -1, 'A': -1, 'I': -1, 'F': -1}
            textList = mainText['documents'][0]['passages']
            abbreDirt = abbreText['documents'][0]['passages']
            abbreList = {}
            for item in abbreDirt:
                abbreList[item['text_short']] = item['text_long_1']
            
            for text in textList:
                section = text['infons']['iao_name_1']
                if section in sectionType:
                    sectionDict[sectionType[section]] += 1
                elif 'iao_name_2' in text['infons']:
                    elsesection = text['infons']['iao_name_2']
                    if elsesection in sectionType:
                        section = elsesection
                        sectionDict[sectionType[section]] += 1
                else:
                    continue
                #print(section)
                offset = int(text['offset'])
                # Find the annotations
                annotation = text['annotations']
                annotList = []
                if len(annotation):
                    for i in annotation:
                        #annotList[i['text']] = i['locations'][0]['offset']
                        annotList.append((i['text'], i['locations'][0]['offset'], i['infons']['identifier']))
                        entity = i['text']
                        identifier = i['text']
                        identifier = identifier.lower().strip(',').strip(';').strip('"').strip("'").strip(':').strip('-').strip(']').strip(')').strip('[').strip('(').strip(',').strip(';').strip('"').strip("'").strip(':').strip('-').rstrip('s')
                        identifier = identifier.strip('')
                        if identifier in enzdict:
                            enzdict[identifier] += 1
                        else:
                            enzdict[identifier] = 1
                        '''
                        if identifier in enzdict:
                            if entity in enzdict[identifier]:
                                enzdict[identifier][entity] += 1
                            else:
                                enzdict[identifier][entity] = 1
                        else:
                            enzdict[identifier] = {entity:1}
                        annonum += 1
                        '''
                
                '''annotList = sorted(annotList, key = lambda x: x[1])
                annonum += len(annotList)
                if len(annotList) == 0:
                    continue
                
                
                # print(annotList)
                # Get the main text
                mainText = text['text']
                doc = nlp(mainText)
                sents = [str(s) for s in doc.sents]
                sentStart = offset
                
                sent_id = -1
                
                for sent in sents:
                    bias = 0     # Replace the abbreviation with defination
                    flag = False # To ensure which sentence has enzymes
                    #print(sent)
                    sent_id += 1
                    sentEnd = sentStart + len(sent)
                    for entity, loc, identifier in annotList:
                        if loc < sentStart:
                            continue
                        # Check if current sentence includes entities
                        elif loc >= sentStart and loc < sentEnd:
                            annotStart = loc - sentStart + bias
                            annotEnd = annotStart + len(entity)
                            # Check if the location is right 
                            if entity == sent[annotStart : annotEnd]:
                                flag = True
                                # Replace the abbreviation using its defination
                                if entity in abbreList:
                                    bias = len(abbreList[entity]) - len(entity)
                                    sent = sent[:annotStart] + sent[annotStart:].replace(entity, abbreList[entity], 1)
                                    entity = abbreList[entity]
                                    annotEnd = annotStart + len(entity)
                                #df_annot = df_annot.append(annotation_output(fileID, sectionType[section], sectionDict[sectionType[section]], sent_id, annotStart, annotEnd, entity), ignore_index = True)
                                
                        else:
                            break
                    #if flag:
                        #df_text = df_text.append(corpus_output(fileID, sectionType[section], sectionDict[sectionType[section]], sent_id, sent.replace('\n', '')), ignore_index = True)
                    #sentStart = sentEnd + 1
                
        df_text.sort_values(['corpus', 'section'], ignore_index = True, inplace = True)
        df_text = df_text.astype({'corpus': str, 'section': str, 'text': str})
        df_text.to_csv(outputText_path, encoding = 'utf_8_sig', sep = '\t', header = False, index = False)
        
        df_annot.sort_values(['corpus', 'section', 'start'], ignore_index=True, inplace=True)
        df_annot = df_annot.astype({'corpus': str, 'section': str, 'start': int, 'end': int, 'enzyme': str})
        df_annot.to_csv(outputAnnot_path, encoding = 'utf_8_sig', sep = '\t', header = False, index = False)
    
    '''
    print("Successfully segment the sentences in {} seconds.".format(time.time() - startTime))
    f = open('D:/BMR-DS/Project 1/Dataset/enzyme_annalysis.json', 'w')
    json.dump(enzdict, f, indent = 4)
    f.close()
       
        
        
        
        
        
    
    
    