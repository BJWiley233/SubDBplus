#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 01:43:17 2020

@author: coyote
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



# query MySQL database to get node and edge infromation
from queries import query_up, query_intact, query_merops2, query_psp2, query_ppase2

userMySQL='dummy'
passwordMySQL='password'
hostMySQL='127.0.0.1'
dbMySQL = 'protTest'

# make cursor dict
cnx = connect(userMySQL, passwordMySQL, hostMySQL)
cursor = cnx.cursor(dictionary=True)
use_database(cnx, cursor, dbMySQL)


# get query for fill neo4j inital UniProt protein load
cursor.execute(query_up)
uniprot_data = cursor.fetchall()
len(uniprot_data)

# get query for fill intact neo4j
cursor.execute(query_intact)
intact_data = cursor.fetchall()
len(intact_data)
#intact_data[987]

# get query for fill merops neo4j
cursor.execute(query_merops2)
merops_data = cursor.fetchall()
len(merops_data)
columns='cleavageID,Ref,cleavage_type,Substrate_formula,\
	substrateName,substrateGenePreferred,\
    subGeneAlt,substrateTax,\
	substrateUniprot,substrateOrganism,\
	protease,code,proteaseUniprot,proteaseGenePreferred,\
    proteaseName,proteaseAltNames,\
    proteaseTaxId,proteaseOrganism,meropsID'.split(",")
merops_df = pd.DataFrame(merops_data, columns=columns)
merops_df['proteaseUniprot']
DUBS = pd.read_csv("/home/coyote/tools/SubDBplus/data/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.03.22-18.34.27.50.tsv",
                   sep="\t")
DUBS['Entry'].values in merops_df['proteaseUniprot'].values
in_it = [i in merops_df['proteaseUniprot'].values for i in DUBS['Entry'].values]
DUBS['Entry'].values[in_it]

# get query for fill psp neo4j
cursor.execute(query_psp2)
psp_data = cursor.fetchall()
len(psp_data)

# get query for fill ppase neo4j
cursor.execute(query_ppase2)
ppase_data = cursor.fetchall()
len(ppase_data)

###############################################################################
# FILL IT!@

# desktop version
#uri = "neo4j://localhost:11003"
# mostly need this host
uri = 'neo4j://localhost:7688'
user = 'neo4j'
# you'll need a password after you start 
password = 'neo4jneo4j'
db = 'protTest'

from fill_neo4j import ProteinExample
prot_db = ProteinExample(uri, user, password)
#prot_db.create_db(db)
#prot_db._drop(db)
prot_db.use_database(db)
prot_db.database

prot_db.create_proteins([uniprot_data[1]])
## test adding meropsID
#prot_db.create_proteins(uniprot_data[0:2])
prot_db.create_intact_interactions(intact_data)
prot_db.create_merops_interactions(merops_data[100:])
prot_db.create_psp_interactions(psp_data)
prot_db.create_depod_interactions(ppase_data)
