#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 20:16:14 2019

@author: sofia
"""
import os
import requests,sys
import argparse
import re
from fuzzywuzzy import process


def main():
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-q', '--question', help='Enter a valid question')
    args = parser.parse_args()
    y = args.question
    a=y.split('<')[0]
    #print(y)
    #print(a)
    

    def choices(x):
        question_list=[
                       (0,'Which gene is upstream to gene <>?'),
                       (0,'upstream to gene <>?'),
                       (1,'How many transcripts does the gene <> have?'), 
                       (1,'How many transcripts of the gene <> ?'),
                       (2,'What is in <>?'),
                       (3,'What is the name of the gene <>?'),
                       (3,'name of the gene <> ?'),
                       (3,'the name of the gene <> ?'),
                       (4,'Where is the location of <>?'),
                       (4,'location of <>?'),
                       (4,'Where is the <>?'),
                       (5,'Which are the refseq transcripts of the gene <>?'), 
                       (5,'refseq transcripts of the gene <>?'),
                       (6,'How many exons does have <>?'), 
                       (6,'How many exons of  <>?'),
                       (7,'Where is the location of the CDS of the <> ?'),
                       (7,'location of the CDS of the <>?'),
                       (8,'What is the primary transcript of the gene <>?'), 
                       (8,'primary transcript of the gene <>?'),
                       (9,'In which gene does belong <>?'),
                       (9,'In which gene is <>?'),
                       (10,'In which pathways is the gene  involved <>?'), 
                       (10,'pathways of the gene involved <>?'),
                       (11,'Which proteins are transcribed from the <> ?'), 
                       (11,'transcribed proteinsof the  <>?'),
                       (12,'What is the function of the gene <> ?'),
                       (12,'function of the gene <>?'),
                       (13,'Which organisms have a homologous gene to <>?'),
                       (14, 'What type is the gene <> ?'), 
                       (14,'type of gene <> ?'),
                       (15, 'What is the sequence in chromosome <> between start(int) and end(int)?'),
                       (15,'sequence in chromosome <> between start(int) and end(int)?'),
                       (16,'What is the allele frequency of <> ?'), 
                       (17,'With which conditions is the mutation associated ?'),
                       (18, 'What is the rs-id of <> ?'), 
                       (18,'rs-id of <> '),
                       (19,'What is the effect of <> ?'),
                       (19,'the effect of <>')]
    
        for question in question_list:
            limit=process.extractOne(x, question_list)
        print(f'Your question is closest to this one: {limit}', limit[0][0])
        return limit[0][0]
    
    result = choices(a)

    if result == 19:
    
        def get_the_entry(y):
            mutation = (re.search(r"\<([A-Za-z0-9_.:>]+)\>", y).group(1))
            return mutation
    
        def get_chrom_transcipt(x): 
    
            url = 'http://mygene.info/v3/query'
            parameters = {
              'fields': 'genomic_pos',
              'species': 'human',
              'q' : 'ensembl.transcript:' + x}
            response = requests.get(url, params=parameters)
            data = response.json()
            return data['hits'][0]['genomic_pos']['chr']
    
        def get_chrom_gene(x): 
            url = 'http://myvariant.info/v1/query'
            parameters = {
                'q': x,
                'fields': 'dbsnp'}
    
            d = (requests.get(url, params=parameters)).json()
            chromosome = d['hits'][0]['dbsnp']['chrom']
            return chromosome
    
    
        def mutation(x):
    
                if 'rs' in x:
                    return x
                elif x.startswith(('chr','Chr', 'cHr','chR', 'CHR', 'cHR', 'ChR')):
                    if '>' in x:
                        return x
                elif x[0].isdigit():
                        return 'chr'+x
                elif 'c.' in x:
                    if x.startswith(('ENST0', 'NM')):
                        position = x.split('.')[1]
                        transcript = re.search(r"([A-Za-z0-9_]+)", choices[0]).group()
                        mutation = 'chr' + get_chrom_transcipt(transcript) + ':' + position
                        return  mutation
                    else:
                        position = x.split('.')[1]
                        gene = re.search(r"([A-Za-z0-9_]+)", choices[0]).group()
                        mutation='chr'+ get_chrom_gene(gene) + ':' + position
                        return  mutation
                else:
                    lista_formats=['rs12345678', 'chr1:1234567A>G', '1:1234567A>G', '<Gene>:c.100A>G', 
                                   '<Transcript (Ensembl or Refeseq_ID)>:c.100A>G']
                    raise Exception("Please specify your mutation properly. The examples for acceptable formats are:" +'\n'+ 
                                    f"{lista_formats}")
    
        k = mutation(get_the_entry(y))
        print(k)
        
        def do_appropriate_mutation_request(k):
    
            if 'rs' in k:
                url = 'http://myvariant.info/v1/query'
                parameters = {
                    'q': k,
                    'fields': 'dbsnp'}
                d = (requests.get(url, params=parameters)).json()
                ID = d['hits'][0]['_id']
                url = 'http://myvariant.info/v1/variant/' + ID
                print(url)
                try: 
                    r = requests.get(url)
                    r.raise_for_status()
    
                except requests.exceptions.HTTPError as errh:
                    print ("Http Error:",errh)
                    #print(f"{data['error']}")
    
                except requests.exceptions.ConnectionError as errc:
                    print ("Error Connecting:",errc)
    
                except requests.exceptions.RequestException as err:
                    print ("OOps: Something Else",err)
    
    
                data = (requests.get(url)).json()
    
                try :
                    effect = data['snpeff']['ann']['effect']
                    print('According to:' + '\n' + f"snpeff: {effect}")
                except KeyError:
                    print(f"There in no snpeff data for effect of this {k} mutation")
    
                try:
                    poly_val = data['cadd']['polyphen']['val']
                    pass
                    poly_cat = data['cadd']['polyphen']['cat']
                    print(f"Polyphen: {poly_val}({poly_cat})")
    
                except KeyError:
                    print(f"There in no cadd data for polyphen of this {k} mutation")
    
                try:
                    cr = data['dbnsfp']['sift']['converted_rankscore']
                    pass    
                    pred = data['dbnsfp']['sift']['pred']
                    score = data['dbnsfp']['sift']['score']
                    print(f"SIFT: converted_rankscore:{cr}, pred:'{pred}', score: {score}")
    
                except KeyError:
                    print(f"There in dbsnf data for this {k} mutation")
    
            else:
    
                url = 'http://myvariant.info/v1/variant/' + k
    
                print(url)
    
                try: 
                    r = requests.get(url)
                    r.raise_for_status()   
    
                except requests.exceptions.HTTPError as errh:
                    print ("Http Error:",errh)
                    #print(f"{data['error']}")
    
                except requests.exceptions.ConnectionError as errc:
                    print ("Error Connecting:",errc)
    
                except requests.exceptions.RequestException as err:
                    print ("OOps: Something Else",err)
    
                    effect = data['snpeff']['ann']['effect']
                    print('According to:' + '\n' + f"snpeff: {effect}")     
                except KeyError:
                    print(f"There in no snpeff data for effect of this {k} mutation")
    
                data = (requests.get(url)).json()
    
                try :
                    effect = data['snpeff']['ann']['effect']
                    print('According to:' + '\n' + f"snpeff: {effect}")
                except KeyError:
                    print(f"There in no snpeff data for effect of this {k} mutation")
    
                try:
                    poly_val = data['cadd']['polyphen']['val']
                    pass
                    poly_cat = data['cadd']['polyphen']['cat']
                    print(f"Polyphen: {poly_val}({poly_cat})")
    
                except KeyError:
                    print(f"There in no cadd data for polyphen of this {k} mutation")
    
                try:
                    cr = data['dbnsfp']['sift']['converted_rankscore']
                    pass    
                    pred = data['dbnsfp']['sift']['pred']
                    score = data['dbnsfp']['sift']['score']
                    print(f"SIFT: converted_rankscore:{cr}, pred:'{pred}', score: {score}")
    
                except KeyError:
                    print(f"There in dbsnf data for this {k} mutation")
    
        do_appropriate_mutation_request(k)
    else:
        print('poulo')

if __name__ == '__main__':
   main()
