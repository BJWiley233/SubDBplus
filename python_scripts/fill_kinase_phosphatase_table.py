#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 20:56:38 2020

@author: coyote
"""

import create_mysql_db
from create_mysql_db import *
import pandas as pd
from io import StringIO
import numpy as np
import json
import re
import requests
import xml.etree.ElementTree as ET


user='dummy'
password='password'
host='127.0.0.1'
db_name = "protTest"

cnx = connect(user, password, host)
cursor = cnx.cursor()
use_database(cnx, cursor, db_name)
cnx.database


TABLES = {}

TABLES['kinasePhosphoSitePlus'] = (
    "CREATE TABLE kinasePhosphoSitePlus ("
    "    uniProtIDKin CHAR(15) NOT NULL,"
    "    geneNamePreferredKin VARCHAR(45),"
    "    protNamePreferredKin VARCHAR(255),"
    "    kinTaxid INT,"
    "    kinOrganism VARCHAR(45),"
    "    uniProtIDSub CHAR(15) NOT NULL,"
    "    geneNamePreferredSub VARCHAR(45),"
    "    geneNameAltSub VARCHAR(45),"
    "    protNamePreferredSub VARCHAR(255),"
    "    subTaxid INT,"
    "    subOrganism VARCHAR(45),"
    "    subModSite CHAR(10),"
    "    sitePlusMinus7AA CHAR(20),"
    "    inVivo BOOLEAN,"
    "    inVitro BOOLEAN,"
    "    UNIQUE(uniProtIDKin, uniProtIDSub, subModSite, inVivo, inVitro)"
    ") ENGINE=InnoDB"
);


TABLES['depodPhosphatase'] = (
    "CREATE TABLE depodPhosphatase ("
    "    uniProtIDPPase CHAR(15) NOT NULL,"
    "    geneNamePreferredPPase VARCHAR(45),"
    "    protNamePreferredPPase VARCHAR(255),"
    "    ppaseTaxid INT,"
    "    ppaseOrganism VARCHAR(45),"
    "    uniProtIDSub CHAR(15) NOT NULL,"
    "    geneNamePreferredSub VARCHAR(45),"
    "    protNamePreferredSub VARCHAR(255),"
    "    subTaxid INT,"
    "    subOrganism VARCHAR(45),"
    "    subDephospoSites VARCHAR(256),"
    "    sitePlusMinus5AA VARCHAR(256),"
    "    inVivo BOOLEAN,"
    "    inVitro BOOLEAN,"
    "    literature VARCHAR(256),"
    "    reliabilityScore INT,"
    "    UNIQUE(uniProtIDPPase, uniProtIDSub, subDephospoSites, inVivo, inVitro)"
    ") ENGINE=InnoDB"
);
#cursor.execute("DROP TABLE IF EXISTS kinasePhosphoSitePlus")
#cursor.execute("DROP TABLE IF EXISTS depodPhosphatase")

for table_name in TABLES:
    table_description = TABLES[table_name]
    try:
        print("Creating table: {}".format(table_name))
        cursor.execute(table_description)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
            print("Table already exists")
        else:
            print(err.msg)


###############################################################################
# PSP
###############################################################################
## Read in PSP data
psp = pd.read_csv("../data/PhosphoSitePlus_Substrates_of_Kinases/Kinase_Substrate_Dataset2", 
                     sep="\t", skiprows=0)
psp.SUB_GENE.fillna(psp.SUBSTRATE, inplace=True)

psp['in_vivo'] = np.where(psp['IN_VIVO_RXN']=='X', True, False)
psp['in_vitro'] = np.where(psp['IN_VITRO_RXN']=='X', True, False)
psp['kin_tax'] = psp['KIN_ORGANISM'].map({'human': 9606,
                                          'rat' : 10116,
                                          'mouse': 10090,
                                          'hamster': 10029})
psp['sub_tax'] = psp['SUB_ORGANISM'].map({'human': 9606,
                                          'rat' : 10116,
                                          'mouse': 10090,
                                          'hamster': 10029})

# only get rows for human, mouse, rat, hamster
psp2 = psp[~psp.filter(like='_tax').isna().any(1)].copy()
psp2['kin_organ'] = psp2['KIN_ORGANISM'].map({'human': 'Homo sapiens',
                                          'rat' : 'Rattus norvegicus',
                                          'mouse': 'Mus musculus',
                                          'hamster': 'Cricetulus griseus'})
psp2['sub_organ'] = psp2['SUB_ORGANISM'].map({'human': 'Homo sapiens',
                                          'rat' : 'Rattus norvegicus',
                                          'mouse': 'Mus musculus',
                                          'hamster': 'Cricetulus griseus'})

psp2['KIN_ACC_ID'] = psp2['KIN_ACC_ID'].str.replace("-.*", "")
psp2['SUB_ACC_ID'] = psp2['SUB_ACC_ID'].str.replace("-.*", "")
psp_final = psp2.replace({np.nan: None})
psp_final["kin_protein_name"] = ""
psp_final["sub_protein_name"] = ""
psp_final[psp_final['SUB_ACC_ID']=='P63058_OBS']
###############################################################################
## Want to get full protein names for the PSP Dataset which didn't come with names
ns = {'up': 'http://uniprot.org/uniprot'}
id_name = dict.fromkeys(pd.concat([psp_final['KIN_ACC_ID'], psp_final['SUB_ACC_ID']]))
len(id_name)
id_name['P63058_OBS']

keys = list(id_name.keys())[0:5]
i = 0
#for key in keys:
    P06537
for key in list(id_name.keys())[i-1:]:
    i += 1
    print("************************", i, key) 
    if id_name[key]:
        continue

    try:
        xml_req = requests.get('https://www.uniprot.org/uniprot/{}.xml'.format(key))
        xml_req.raise_for_status()
        print("good")
        if xml_req.content:
            tree = ET.fromstring(xml_req.content) 
            if tree.find("./up:entry/up:protein/up:recommendedName/up:fullName", ns) is not None:
                print(1, key)
                fullName = tree.find("./up:entry/up:protein/up:recommendedName/up:fullName", ns).text
            elif tree.find("./up:entry/up:protein/up:submittedName/up:fullName", ns) is not None:
                print(2, key)
                fullName = tree.find("./up:entry/up:protein/up:submittedName/up:fullName", ns).text
            else:
                print(3, key)
                fullName = None
            id_name[key] = fullName
        else:
            print("Bad request for", key)
       
    except requests.exceptions.HTTPError as e:
        print("Error getting XML for {} ".format(key))
        print(e)
        

psp_final["kin_protein_name"] = psp_final.apply(lambda x: id_name[x['KIN_ACC_ID']], axis=1)
psp_final["sub_protein_name"] = psp_final.apply(lambda x: id_name[x['SUB_ACC_ID']], axis=1)

###############################################################################
#Enter kinases in MySQL
kinase_log = "kinase_not_entered.log"
open(kinase_log, 'w').close()

## SUB_GENE (column 8 in data) = geneNamePreferredSub: 
    ## the lowercase for mouse/rat (good!)
## SUBSTRATE (column 5 in data) = geneNameAltSub: 
    ## alternate name or sometimes uppercase for mouse/rat (don't like uppercase if mouse/rat)

for idx, row in psp_final.iterrows():
    print(idx)
    try:
        # uniProtIDKin, geneNamePreferredKin, protNamePreferredKin, kinTaxid, kinOrganism
        # uniProtIDSub, geneNamePreferredSub, geneNameAltSub
        cursor.execute(
            "REPLACE INTO kinasePhosphoSitePlus VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
            (row['KIN_ACC_ID'], row['GENE'], row['kin_protein_name'], row['kin_tax'], row['kin_organ'], 
             row['SUB_ACC_ID'], row['SUB_GENE'], row['SUBSTRATE'], row['sub_protein_name'],
             row['sub_tax'], row['sub_organ'], row['SUB_MOD_RSD'],
             row['SITE_+/-7_AA'], row['in_vivo'], row['in_vitro'])
        )
        cnx.commit()
    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))
        with open(kinase_log, 'a') as f:
            f.write(f"{idx+1}:    {row['KIN_ACC_ID']}    {row['SUB_ACC_ID']}\n")





###############################################################################
# DEPOD
###############################################################################
## Read in DEPOD data
ppase = pd.read_csv("../data/PPase_protSubtrates_201903.csv", na_values=["N/A", "N/A "])
ppase.columns
#ppase = ppase.replace({np.nan: None})

ppase['in_vivo'] = np.where(ppase['Assay/s to infer dephosphorylation'].str.match(".*vivo"), True, False)
ppase['in_vitro'] = np.where(ppase['Assay/s to infer dephosphorylation'].str.match(".*vitro"), True, False)
ppase[['Assay/s to infer dephosphorylation', 'in_vivo','in_vitro']]
ppase_final = ppase.replace({np.nan: None})
######################
# how many sites known?
sum([len(i) for i in ppase_final['Dephosphosites'].str.split(",") if i is not None])
######################
ppase_final["phos_protein_name"] = ""
ppase_final["sub_protein_name"] = ""

###############################################################################
## Want to get full protein names for the DEPOD Dataset which didn't come with names
ns = {'up': 'http://uniprot.org/uniprot'}
id_name_ppase = dict.fromkeys(pd.concat([ppase_final['Phosphatase accession numbers'], 
                                         ppase_final['Substrate accession numbers']]))
len(id_name_ppase)

i = 0
#for key in keys:
for key in id_name_ppase.keys():
    i += 1
    print("************************", i, key) 
    
    try:
        xml_req = requests.get('https://www.uniprot.org/uniprot/{}.xml'.format(key))
        xml_req.raise_for_status()
        print("good")
        tree = ET.fromstring(xml_req.content) 
        if tree.find("./up:entry/up:protein/up:recommendedName/up:fullName", ns) is not None:
            print(1, key)
            fullName = tree.find("./up:entry/up:protein/up:recommendedName/up:fullName", ns).text
        elif tree.find("./up:entry/up:protein/up:submittedName/up:fullName", ns) is not None:
            print(2, key)
            fullName = tree.find("./up:entry/up:protein/up:submittedName/up:fullName", ns).text
        else:
            print(3, key)
            fullName = None
        id_name_ppase[key] = fullName
    except requests.exceptions.HTTPError as e:
        print("Error getting XML for {} ".format(key))
        print(e)

ppase_final["phos_protein_name"] = ppase_final.apply(lambda x: id_name_ppase[x['Phosphatase accession numbers']], axis=1)
ppase_final["sub_protein_name"] = ppase_final.apply(lambda x: id_name_ppase[x['Substrate accession numbers']], axis=1)

for idx, row in ppase_final.iterrows():
    ppase_final['phos_protein_name'] = id_name_ppase[row['Phosphatase accession numbers']]
    row['sub_protein_name'] = id_name_ppase[row['Substrate accession numbers']]
    
###############################################################################
#Enter phosphatases in MySQL
phosphatase_log = "phos_not_entered.log"
open(phosphatase_log, 'w').close()

for idx, row in ppase_final.iterrows():
    print(idx)
    try:
        cursor.execute(
            "REPLACE INTO depodPhosphatase VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
            (row['Phosphatase accession numbers'], row['Phosphatase entry names'], 
             row['phos_protein_name'], 9606, 'Homo sapiens', row['Substrate accession numbers'], 
             row['Substrate entry names'], row['sub_protein_name'], 9606, 'Homo sapiens',
             row['Dephosphosites'], row['5 amino acid window around the dephosphosite (small letters)'],
             row['in_vivo'], row['in_vivo'], row['Literature source/s'], row['Reliability score'])
        )
        cnx.commit()
    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))
        with open(phosphotase_log, 'a') as f:
            f.write(f"{idx+1}:    {row['Phosphatase accession numbers']}    {row['Substrate accession numbers']}\n")
    
