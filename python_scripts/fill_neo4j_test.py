#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 05:10:10 2020

@author: coyote

This is to test overlapping nodes in TntaAt with the other databases
"""

import logging
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import create_mysql_db
from create_mysql_db import *
import pandas as pd
import numpy as np
import requests
import re
import json

############################################################################
"""
This script also loads data into the MySQL data warehouse
based on the issue that there was not gene name data for
the MEROPS substrates for the 'name' or label of the nodes
in Neo4j
"""
############################################################################


userMySQL='root'
passwordMySQL='**'
hostMySQL='127.0.0.1'
dbMySQL = 'protTest'

cnx = connect(userMySQL, passwordMySQL, hostMySQL)
cursor = cnx.cursor(dictionary=True)
use_database(cnx, cursor, dbMySQL)

## testing data and queries for loading into neo4j
from queries import query_merops
cursor.execute(query_merops)
merops_sel = cursor.fetchall()
merops = pd.DataFrame([[i['interactorA'], i['interactorB'], 
                        i["interactionID"], i['uniProtID'],
                        i['cleavageID']] for i in merops_sel])
merops
merops.columns = ['interactorA','interactorB',"interactionID", 'uniProtID', 'cleavageID']
merops_IDs = np.unique(merops['cleavageID'])

from queries import query_psp
cursor.execute(query_psp)
psp_sel = cursor.fetchall()
psp = pd.DataFrame([[i['interactorA'], i['interactorB'], 
                     i["interactionID"], i['uniProtID']] for i in psp_sel])
psp
psp.columns = ['interactorA','interactorB',"interactionID", 'uniProtID']


from queries import query_ppase
cursor.execute(query_ppase)
ppase_sel = cursor.fetchall()
ppase = pd.DataFrame([[i['interactorA'], i['interactorB'], 
                       i["interactionID"], i['uniProtID']] for i in ppase_sel])
ppase
ppase.columns = ['interactorA','interactorB',"interactionID", 'uniProtID']

np.unique([merops['interactorA'], merops['interactorB']])
np.unique([psp['interactorA'], psp['interactorB']])
np.unique(psp['interactorB'])
np.unique([ppase['interactorA'], ppase['interactorB']])
np.unique(ppase['interactorB'])

get_uniprot = np.unique((merops['interactorA'].tolist() + 
                         merops['interactorB'].tolist() +
                         psp['interactorA'].tolist() +
                         psp['interactorB'].tolist() +
                         ppase['interactorA'].tolist() + 
                         ppase['interactorB'].tolist() +
                         merops['uniProtID'].tolist() +
                         psp['uniProtID'].tolist() +
                         ppase['uniProtID'].tolist()))

get_intact = np.unique((merops['interactionID'].tolist() + 
                         psp['interactionID'].tolist() +
                         ppase['interactionID'].tolist()))

len(get_uniprot)
len(get_intact)

###############################################################################
## Not getting good data for gene names for merops substrates
## Need to get merops substrate genes names from UniProt
## Merops only had a description column which was hard to 
## get gene names for nodes. Do this with batch request

query_merops3 = """
SELECT DISTINCT s.Uniprot AS substrateUniprot
FROM protTest.Substrate_search s
INNER JOIN 
	protTest.proteinsUniprot p
    ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%%', p.organism, '%%')
WHERE s.Uniprot IS NOT NULL
"""
## non dict cursor
cursor = cnx.cursor()
cursor.execute(query_merops3)
merops_subs = cursor.fetchall()
arr = [i[0] for i in merops_subs]
len(merops_subs)
merops_subs2 = [i[0][0:6] for i in merops_subs]
up_tbl_cols = ["id", "reviewed", "entry name", "genes(PREFERRED)",
               "genes(ALTERNATIVE)", "organism-id", "organism"]


# Uniprot batch request docs!
#https://www.biostars.org/p/304247/
#https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf
BASE = "http://www.uniprot.org"
TOOL_ENDPOINT = "/uploadlists/"
def map_retrieve(ids2map, source_fmt='ACC+ID',target_fmt='ACC', output_fmt='tab'):
    if hasattr(list(ids2map), 'pop'):
        ids2map = ' '.join(ids2map)
 
    payload = {'query': ids2map,
               'format': output_fmt,
               'columns': ','.join(up_tbl_cols),
               'from': source_fmt,
               'to': target_fmt}
    
    response = requests.get(BASE + TOOL_ENDPOINT, params=payload)
    if response.ok:
        protein_df = pd.read_csv(StringIO(response.text), sep="\t")
        ## pandas is good but need to pass None to MySQL
        df_final = protein_df.replace({np.nan: None})
        df_final2 = df_final.iloc[:, :-1]
        df_final2.columns
        return df_final2
    else:
        print(response.raise_for_status())

chunks = np.array_split(merops_subs2, 10)
frames = [map_retrieve(i) for i in chunks]
merops_up_sub_tbl = pd.concat(frames)
reviewed = merops_up_sub_tbl[merops_up_sub_tbl['Status'] == 'reviewed']

## add new table to database for merops substrate gene names
TABLES = {}
TABLES['meropsSubsUniprot'] = """
CREATE TABLE meropsSubsUniprot (
    uniProtID CHAR(15) NOT NULL PRIMARY KEY,
    entryName VARCHAR(45) NOT NULL,
    geneNamePreferred VARCHAR(45),
    geneNamesAlternative JSON DEFAULT NULL,
    taxid INT,
    organism VARCHAR(100)
    ) ENGINE=InnoDB
"""
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

log_file = "uniprot_merops_subs_not_entered.log"
open(log_file, 'w').close()

for idx, row in reviewed.iterrows():
#for idx, row in row.iterrows():    
    print("*****************************",  row['Entry'])
    uniProtID = row['Entry']
    entryName = row['Entry name']
    
    if row['Gene names  (primary )']:
        geneNamePreferred = row['Gene names  (primary )'].split(";")[0]
    else:
        geneNamePreferred = None
   
    ## can be None
    if row['Gene names  (synonym )']==None:
        geneNamesAlternative=None
    else:
        geneNamesAlternative = row['Gene names  (synonym )'].split(' ')
        
    organismID = row['Organism ID']
    organism = row['Organism'].split(" (")[0]
    geneNamesAlternative_ = None if geneNamesAlternative==None else json.dumps(geneNamesAlternative)

    
    try:
        cursor.execute(
            "INSERT IGNORE INTO meropsSubsUniprot VALUES (%s,%s,%s,%s,%s,%s)",
            (uniProtID, entryName, geneNamePreferred, geneNamesAlternative_,  
             organismID, organism))
        cnx.commit()
    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))
        with open(log_file, 'a') as f:
            f.write(f"{row['Entry']}\n")


###############################################################################

# Test filling neo4j

# mostly need this host
uri = 'neo4j://localhost:7687'
user = 'neo4j'
# you'll need a password after you start 
password = '**'
db = 'protTest'


from fill_neo4j import ProteinExample
prot_db = ProteinExample(uri, user, password)
#prot_db.create_db(db)
#prot_db._drop(db)
prot_db.use_database(db)
prot_db.database

# query MySQL database to get node and edge infromation
# Need to remove comments in queries.py if testing again
from queries import query_up, query_intact, query_merops2, query_psp2, query_ppase2

#https://stackoverflow.com/questions/589284/imploding-a-list-for-use-in-a-python-mysqldb-in-clause
format_strings = ','.join(['%s'] * len(get_uniprot))
cursor.execute(query_up % format_strings, tuple(get_uniprot))
uniprot_data = cursor.fetchall()
len(uniprot_data)

# get query for fill intact neo4j
format_strings_intact = ','.join(['%s'] * len(get_intact))
cursor.execute(query_intact % format_strings_intact, tuple(get_intact))
intact_data = cursor.fetchall()
len(intact_data)

# get query for fill merops neo4j
format_strings_merops = ','.join(['%s'] * len(merops_IDs))
cursor.execute(query_merops2 % format_strings_merops, tuple(merops_IDs))
merops_data = cursor.fetchall()
len(merops_data)

# get query for fill psp neo4j
format_strings_psp = ','.join(['%s'] * len(np.unique(psp['interactorB'])))
cursor.execute(query_psp2 % format_strings_psp, tuple(np.unique(psp['interactorB'])))
psp_data = cursor.fetchall()
len(psp_data)

# get query for fill ppase neo4j
format_strings_ppase = ','.join(['%s'] * len(np.unique(ppase['interactorB'])))
cursor.execute(query_ppase2 % format_strings_ppase, tuple(np.unique(ppase['interactorB'])))
ppase_data = cursor.fetchall()
len(ppase_data)


"""
Tests
prot_db.create_proteins(uniprot_data)
prot_db.create_intact_interactions(intact_data)
prot_db.create_merops_interactions(merops_data)
prot_db.create_psp_interactions(psp_data)
prot_db.create_depod_interactions(ppase_data)
"""
