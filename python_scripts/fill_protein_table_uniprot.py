#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 04:53:27 2020

@author: coyote
"""


import create_mysql_db
from create_mysql_db import *

import pandas as pd
from io import StringIO
import requests
import numpy as np
import json
import xml.etree.ElementTree as ET
import re



# humans, mouse, rat, chinese hamster, Escherichia coli, Saccharomyces cerevisiae, 
# Arabidopsis, drosophila
# 
organisms = [9606, 10090, 10116, 10029, 83333, 559292, 3702, 7227]
# growth factors??
families = ['kinase', 'phosphatase', 'transferase', 'methyltransferase', 
            'peptidase', 'isomerase', 'G-protein coupled receptor']
# https://www.uniprot.org/help/uniprotkb_column_names
# Real tables names on website and download are different
# example "Gene names (primary)" for "genes(PREFERRED)"
up_tbl_cols = ["id", "entry name", "reviewed", "protein names", "genes(PREFERRED)",
               "genes(ALTERNATIVE)", "genes" , "interactor", "organism-id", "organism", 
               "families", "length", "database(MEROPS)"]


user='root'
password='**'
host='127.0.0.1'
db_name = "protTest"

cnx = connect(user, password, host)
cursor = cnx.cursor()

create_database(cursor, db_name)
use_database(cnx, cursor, db_name)

#cursor.execute("DROP TABLE IF EXISTS proteins")

TABLES = {}
TABLES['proteinsUniprot'] = (
    "CREATE TABLE proteinsUniprot ("
    "	 uniProtID CHAR(15) NOT NULL PRIMARY KEY,"
    "    entryName VARCHAR(45) NOT NULL,"
    "    reviewed BOOL NOT NULL,"
    "    proteinName MEDIUMTEXT NOT NULL,"
    "    alternateNames JSON DEFAULT NULL,"
    "    geneNamePreferred VARCHAR(45),"
    "    geneNamesAlternative JSON DEFAULT NULL,"
    "    geneNamesAll JSON DEFAULT NULL,"
    "    interactsWith JSON DEFAULT NULL,"
    "    taxid INT,"
    "    organism VARCHAR(100),"
    "    organismCommon MEDIUMTEXT,"
    "    proteinFamilies VARCHAR(255),"
    "    length INT,"
    "    meropsID VARCHAR(8) DEFAULT NULL,"
    "    headProteinFamily VARCHAR(45) NOT NULL"
    ") ENGINE=InnoDB"
);

TABLES['pdb'] = (
    "CREATE TABLE pdb ("
    "	 pdbID CHAR(5) NOT NULL PRIMARY KEY,"
    "    method VARCHAR(45) NOT NULL"
    ") ENGINE=InnoDB"
);

TABLES['uniprotPdbJoin'] = (
    "CREATE TABLE uniprotPdbJoin ("
    "	 uniProtID CHAR(15) NOT NULL,"
    "    pdbID CHAR(5) NOT NULL,"
    "    chain VARCHAR(25) NOT NULL,"
    "    UNIQUE (uniProtID,pdbID,chain),"
	"	 FOREIGN KEY (uniProtID)" 
    "       REFERENCES proteinsUniprot(uniProtID),"
	"	 FOREIGN KEY (pdbID)"
    "       REFERENCES pdb(pdbID)"
    ") ENGINE=InnoDB"
);

TABLES['meropsSubsUniprot']
for table in TABLES:
    drop_table(cursor, db_name, table)


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
            
    
# https://www.uniprot.org/help/api%5Fqueries
uniprot_url = "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+family:{}+AND+("
uniprot_url += "+OR+".join(["organism:{}".format(str(i)) for i in organisms])
uniprot_url += ")&columns=" + ",".join(up_tbl_cols) + "&format=tab"
## for XML namespace for alt names and pdbs
ns = {'up': 'http://uniprot.org/uniprot'}

urls = [uniprot_url.format(fam) for fam in families]

log_file = "uniprot_not_entered.log"
pdb_log = "pdb_not_entered.log"
open(log_file, 'w').close()
open(pdb_log, 'w').close()


cursor.execute("SELECT DISTINCT uniProtID FROM proteinsUniprot")
uniprot_data = cursor.fetchall()
in_up_tble = [i[0] for i in uniprot_data]

for i in range(0, len(urls)):
#for i in range(0, 1):    
    req = requests.get(urls[i])
    if req.ok:
        protein_df = pd.read_csv(StringIO(req.text), sep="\t")
        ## pandas is good but need to pass None to MySQL
        df_final = protein_df.replace({np.nan: None})
    else:
        req.raise_for_status()
        continue
    

    n=0
    
    for idx, row in df_final.iterrows():
        print("*****************************",  row['Entry'])
        # if already in table skip it
        uniProtID = row['Entry']
        if uniProtID in in_up_tble:
            print(uniProtID, "in table already")
            continue
        
        entryName = row['Entry name']
        reviewed = 1 if row['Status']=='reviewed' else 0
        
        ## no good way to get protein names separate from alternative names
        ## in the tabular format as some names have `tRNA (something) rest of first name`
        ##
        ## if request for XML fails, set protein proteinName to long string
        ## that includes the alternateNames and set alternateNames=NONE
        pdbs = {}
        pdb_chain = {}
        try:
            xml_req = requests.get('https://www.uniprot.org/uniprot/{}.xml'.format(uniProtID))
            xml_req.raise_for_status()
            print("good")
            tree = ET.fromstring(xml_req.content)
            proteinName = tree.find("./up:entry/up:protein/up:recommendedName/up:fullName", ns).text
            alternateNames = [i.text for i in tree.findall("./up:entry/up:protein/up:alternativeName/up:fullName", ns)] 
            if len(alternateNames) > 0:
                pass
            else: 
                alternateNames = None
            ######################################################
            pdb_structs = tree.findall("./up:entry/up:dbReference[@type='PDB']", ns)
            for p in pdb_structs:
                id_ = p.get('id')
                for prop in p.findall("./up:property", ns):
                    if prop.get('type') == 'method':
                        pdbs[id_] = prop.get('value')
                    elif prop.get('type') == 'chains':
                        chains = re.sub("=.*", "", prop.get('value'))
                        pdb_chain[id_] = chains       
            ######################################################
            
        ## default to protein names 
        except requests.exceptions.HTTPError as e:
            print("Error getting XML for {} \
                  defaulting to long entry for ".format(uniProtID))
            print(e)
            proteinName = row['Protein names']
            alternateNames = None
        
        ## can be None
        geneNamePreferred = row['Gene names  (primary )']
       
        ## can be None
        if row['Gene names  (synonym )']==None:
            geneNamesAlternative=None
        else:
            geneNamesAlternative = row['Gene names  (synonym )'].split(' ')
        
        ## won't be None???? can be None
        if row['Gene names']==None:
            geneNamesAll=None
        else:
            geneNamesAll = row['Gene names'].split(' ')
        
        ## can be None
        if row['Interacts with']==None:
            interactsWith=None
        else:
            interactsWith = row['Interacts with'].split('; ')
        
        organismID = row['Organism ID']
        organism = row['Organism'].split(" (")[0]
        organismCommon = "; ".join(re.findall('\(([^)]+)', row['Organism']))
        proteinFamilies = row['Protein families']
        length = row['Length']
        meropsID = None if row['Cross-reference (MEROPS)'] == None else row['Cross-reference (MEROPS)'].split(";")[0]
        headProteinFamily = families[i]
        
        ## slight chance primary gene name is empty
        if not geneNamePreferred and geneNamesAll:
            geneNamePreferred = geneNamesAll[0]
        
        ## picky: json_dumps of None is 'null' not NULL
        alternateNames_ = None if alternateNames==None else json.dumps(alternateNames)
        geneNamesAlternative_ = None if geneNamesAlternative==None else json.dumps(geneNamesAlternative)
        geneNamesAll_ = None if geneNamesAll==None else json.dumps(geneNamesAll)
        interactsWith_ = None if interactsWith==None else json.dumps(interactsWith)
        
        try:
            cursor.execute(
                "INSERT IGNORE INTO proteinsUniprot VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
                (uniProtID, entryName, reviewed, proteinName, alternateNames_, geneNamePreferred,
                 geneNamesAlternative_, geneNamesAll_, interactsWith_, 
                 organismID, organism, organismCommon, proteinFamilies, length, 
                 meropsID, headProteinFamily))
            cnx.commit()
        except mysql.connector.Error as err:
            print("Something went wrong in proteinsUniprot table: {}".format(err))
            with open(log_file, 'a') as f:
                f.write(f"{row['Entry']}\n")
        
        if len(pdbs) > 0:
            for k, v in pdbs.items():
                    try:
                        cursor.execute(
                            "INSERT IGNORE INTO pdb VALUES (%s,%s)", (k, v))
                        cnx.commit()
                    except mysql.connector.Error as err:
                        print("Something went wrong in pdb table: {}".format(err))
                        print(k, v)
                        with open(pdb_log, 'a') as f:
                            f.write(f"{k}\n")
                        continue
        if len(pdb_chain) > 0:
            for k, v in pdb_chain.items():
                try:
                    cursor.execute(
                        "INSERT IGNORE INTO uniprotPdbJoin VALUES (%s,%s,%s)", 
                        (uniProtID, k, v))
                    cnx.commit()
                except mysql.connector.Error as err:
                    print("Something went wrong in uniprotPdbJoin table: {}".format(err))
                    print(uniProtID, k, v)
                    with open(pdb_log, 'a') as f:
                        f.write(f"{k}:   {v}\n")  
                    continue
        n += 1 
        #if n == 500:
            #break

    










