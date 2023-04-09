#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 08:10:50 2020

@author: coyote
"""

import create_mysql_db
from create_mysql_db import *
import pandas as pd
from io import StringIO
import numpy as np
import json
import re



"""
preprocess with R
[1] "X.ID.s..interactor.A" "ID.s..interactor.B"   "geneNamePreferredA"   "geneNamePreferredB"  
[5] "altNamesA"            "altNamesB"            "taxidA"               "taxidB"              
[9] "organismA"            "organismB"            "detectionMethod"      "interactionType"     
[13] "interactionID"        "biologicalRoleA"      "biologicalRoleB"      "pubmedID"            
[17] "sourceDB"             "isNegative"           "direction"  


"""

user='dummy'
password='password'
host='127.0.0.1'
db_name = "protTest"

cnx = connect(user, password, host)
cursor = cnx.cursor()

#create_database(cursor, db_name)
use_database(cnx, cursor, db_name)
cnx.database

cursor.execute("DROP TABLE IF EXISTS intact")
TABLES = {}
## some geneNamePreferred are long because DNA ex.
## psi-mi:cacaccacgttgactaccgt_minus_1-4_nucleotides_at_3'_end(display_short)
## gives cacaccacgttgactaccgt_minus_1-4_nucleotides_at_3'_end
TABLES['intact'] = (
    "CREATE TABLE intact  ("
    "    interactionID CHAR(15) NOT NULL,"
    "	 interactorA VARCHAR(25) NOT NULL,"
    "    interactorB VARCHAR(25),"
    "    geneNamePreferredA MEDIUMTEXT,"
    "    geneNamePreferredB MEDIUMTEXT,"
    "    alternateNamesA JSON DEFAULT NULL,"
    "    alternateNamesB JSON DEFAULT NULL,"
    "    geneNamesAlternativeA JSON DEFAULT NULL,"
    "    geneNamesAlternativeB JSON DEFAULT NULL,"
    "    taxidA INT,"
    "    taxidB INT,"
    "    organismA VARCHAR(100),"
    "    organismB VARCHAR(100),"
    "    detectionMethod TINYTEXT,"
    "    interactionType TINYTEXT NOT NULL,"
    "    biologicalRoleA VARCHAR(45) DEFAULT NULL,"
    "    biologicalRoleB VARCHAR(45) DEFAULT NULL,"
    "    publicationID VARCHAR(255) DEFAULT NULL,"
    "    sourceDB VARCHAR(45) DEFAULT NULL,"
    "    isNegative BOOLEAN DEFAULT FALSE,"
    "    direction CHAR(5) NOT NULL"
    ") ENGINE=InnoDB"
);

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



intact = pd.read_csv("../data/intact_cleaned_R.txt", sep="\t", na_values="NA")

df_final = intact.replace({np.nan: None})
len(df_final)
df_final_unique = df_final.drop_duplicates()
len(df_final_unique)
## Fixed CHEBI for B and A in R script but only for B in MySQL and Neo4j
## B meaning the column from IntAct not in the reaction direction
A = [itm[0] for itm in df_final_unique['X.ID.s..interactor.A'].str.findall('CHEBI.*') if len(itm)>0]
len(A)
df_final_unique[df_final_unique['ID.s..interactor.B'].str.match('CHEBI.*')==True]
##################################################################
len(df_final_unique)
## Found out these CHEBI molecules could be good to link to CHEBI
## 12/3/2020 Fixed the R data cleaning script
## These all zero now
A = [itm[0] for itm in df_final_unique['X.ID.s..interactor.A'].str.findall('^[0-9]{5,6}$') if len(itm)>0]
len(A)
len(np.unique(A))
hasB = df_final_unique['ID.s..interactor.B'][df_final_unique['ID.s..interactor.B'].notnull()]
B = [itm[0] for itm in hasB.str.findall('^[0-9]{5,6}$') if len(itm)>0]
len(B)
## test updating these in Neo4j with Cypher csv query
len(np.unique(B))
len(np.unique(np.concatenate([A, B])))

##################################################################
## Enter into MySQL
log_file = "intact_not_entered.log"
open(log_file, 'w').close()
n=0
for idx, row in df_final_unique.iterrows():
    altNamesA = None if not row["altNamesA"] else row["altNamesA"].split('|')
    altNamesB = None if not row["altNamesB"] else row["altNamesB"].split('|')
    
    alternateNamesA = None if not altNamesA else list(filter(re.compile(" ").search, altNamesA))
    alternateNamesA = None if (not alternateNamesA or len(alternateNamesA)==0) else alternateNamesA
    geneNamesAlternativeA = None if not altNamesA else list(filter(lambda v: not re.compile(" ").search(v), altNamesA))
    geneNamesAlternativeA = None if (not geneNamesAlternativeA or len(geneNamesAlternativeA)==0) else geneNamesAlternativeA
    s
    alternateNamesB = None if not altNamesB else list(filter(re.compile(" ").search, altNamesB))
    alternateNamesB = None if (not alternateNamesB or len(alternateNamesB)==0) else alternateNamesB
    geneNamesAlternativeB = None if not altNamesB else list(filter(lambda v: not re.compile(" ").search(v), altNamesB))
    geneNamesAlternativeB = None if (not geneNamesAlternativeB or len(geneNamesAlternativeB)==0) else geneNamesAlternativeB
    
    ## picky: json_dumps of None is 'null' not NULL
    alternateNamesA_ = None if alternateNamesA==None else json.dumps(alternateNamesA)
    geneNamesAlternativeA_ = None if geneNamesAlternativeA==None else json.dumps(geneNamesAlternativeA)
    alternateNamesB_ = None if alternateNamesB==None else json.dumps(alternateNamesB)
    geneNamesAlternativeB_ = None if geneNamesAlternativeB==None else json.dumps(geneNamesAlternativeB)
    
    if row["interactionType"] == "phosphorylation reaction":
        interactionType = "phosphorylation"
    elif row["interactionType"] == "dephosphorylation reaction":
        interactionType = "dephosphorylation"
    elif row["interactionType"] == "protein cleavage":
        interactionType = "proteolysis"
    ## if cleavage and uniprot ids for both 
    elif row["interactionType"] == "cleavage reaction":
        if len(row["X.ID.s..interactor.A"])==6 and len(row["X.ID.s..interactor.A"])==6:
            interactionType = "proteolysis"     
    else:
        interactionType = row["interactionType"]
    
    try:
        cursor.execute(
            "REPLACE INTO intact VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
            (row["interactionID"], row["X.ID.s..interactor.A"], row["ID.s..interactor.B"],
             row["geneNamePreferredA"], row["geneNamePreferredB"], alternateNamesA_,
             alternateNamesB_, geneNamesAlternativeA_, geneNamesAlternativeB_,
             row["taxidA"], row["taxidB"], row["organismA"], row["organismB"],
             row["detectionMethod"], interactionType, row["biologicalRoleA"], 
             row["biologicalRoleB"], row["pubmedID"], row["sourceDB"], 
             row["isNegative"], row["direction"]))
        cnx.commit()
    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))
        with open(log_file, 'a') as f:
            f.write(f"{row['interactionID']}\n")
    
    print(n)
    n += 1


