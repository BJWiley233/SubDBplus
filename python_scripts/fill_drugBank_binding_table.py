#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 10:18:55 2020

@author: coyote
"""

import create_mysql_db
from create_mysql_db import *
from io import StringIO
import requests
import json
import xml.etree.ElementTree as ET
import re


user='dummy'
#password='s3kr1t' # standard fake password
password='password'
host='127.0.0.1'
db_name = "protTest"

cnx = connect(user, password, host)
cursor = cnx.cursor()
use_database(cnx, cursor, db_name)

#cursor.execute("DROP TABLE IF EXISTS pdbBindingSites")

TABLES = {}
TABLES['pdbDrugBank'] = (
    "CREATE TABLE pdbDrugBank ("
    "    pdbID CHAR(5) NOT NULL,"
    "	 drugBankID CHAR(15) NOT NULL,"
    "    drugShort CHAR(5),"
    "    drugLong VARCHAR(256) NOT NULL,"
    "    UNIQUE(pdbID,drugBankID,drugShort,drugLong),"
    "    FOREIGN KEY (pdbID)"
    "       REFERENCES pdb(pdbID)"
    ") ENGINE=InnoDB"
);

TABLES['pdbBindingSites'] = (
    "CREATE TABLE pdbBindingSites ("
    "    pdbID CHAR(5) NOT NULL,"
    "    siteID CHAR(5),"
    "    structResidNum INT NOT NULL,"
    "    uniprotResidNum INT NOT NULL,"
    "    residue CHAR(12) NOT NULL,"
    "    residChain CHAR(2) NOT NULL,"
    "    ligandResidNum INT,"
    "    ligandShort CHAR(5),"
    "    ligandLong TINYTEXT,"
    "    ligandChain CHAR(2),"
    "    UNIQUE(pdbID,siteID,structResidNum,uniprotResidNum,residue,residChain,ligandResidNum),"
    "    FOREIGN KEY (pdbID)"
    "       REFERENCES pdb(pdbID)"
    ") ENGINE=InnoDB"
);


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
            

## for XML namespace for alt names and pdbs
ns = {'pc': 'http://www.ncbi.nlm.nih.gov'}


drug_log_file = "drug_not_entered.log"
binding_site_log = "bs_not_entered.log"
open(drug_log_file, 'w').close()
open(binding_site_log, 'w').close()

## get all PDBs
cursor.execute("SELECT pdbID FROM pdb;")
pdbids = cursor.fetchall()
## test
#pdbids = ['1A16']

n=0
# drugbank molecules associated with structure
for pdb in pdbids:
    print(pdb['pdbID'])
    pdb_dbUrl = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/drugbank/%s" % pdb['pdbID'];
    req = requests.get(pdb_dbUrl)
    if req.ok:
        print("OK")
    else:
        with open(drug_log_file, 'a') as f:
            f.write(f"{pdb}:\n")
        n += 1
        continue
    
    drug_json = json.loads(req.content)
    for key in drug_json:
        for drugs in drug_json[key]:
            for drug in drugs:
                drugbank_id = drugs[drug]['drugbank_id']
                pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/drugbank/%s/XML" %  drugbank_id
                try:
                    pubchem_req = requests.get(pubchem_url)
                    pubchem_req.raise_for_status()
                    print("good")
                    tree = ET.fromstring(pubchem_req.content)
                    synonym = tree.find(".//pc:PC-Substance_synonyms_E", ns).text
                    
                    try:
                        cursor.execute(
                            "REPLACE INTO pdbDrugBank VALUES (%s,%s,%s,%s)",
                            (pdb['pdbID'], drugbank_id, drug, synonym)
                        )
                        cnx.commit()
                    except mysql.connector.Error as err:
                        print("Something went wrong: {}".format(err))
                        with open(drug_log_file, 'a') as f:
                            f.write(f"{pdb}:  {drugbank_id}\n")
                except requests.exceptions.HTTPError as e:
                    print("Error getting XML for %s" % drugbank_id)
                    print(e)
    
    n += 1 
    print(n)
    if n == 2:
        break
    
n=300
# binding sites in the .pdb files of structures
for pdb in pdbids[n:]:
    
    pdb_ligandUrl = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/%s" % pdb;
    req_lig = requests.get(pdb_ligandUrl)
    if req_lig.ok:
        print("OK")
    else:
        with open(binding_site_log, 'a') as f:
            f.write(f"{pdb}:\n")
        continue
    
    ligands = {}
    lig_json = json.loads(req_lig.content)
    for key in lig_json:
        for lig in lig_json[key]:
            ligands[str(lig['author_residue_number'])] = {
                    'shortName': lig['chem_comp_id'],
                    'longName': lig['chem_comp_name'],
                    'chainID': lig['chain_id']
                }    
    
    pdb_bsUrl = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/%s" % pdb;
    req_bs = requests.get(pdb_bsUrl)
    if req_bs.ok:
        print("OK")
    else:
        with open(binding_site_log, 'a') as f:
            f.write(f"{pdb}:\n")
        continue
    
    bind_json = json.loads(req_bs.content)
    for key in bind_json:
        for site in bind_json[key]:
            site_ID = site['site_id']
            # REMARK 800 SITE_IDENTIFIER: NUL
            if site_ID.upper() == "NUL":
                site_ID = None
            ligand_res_num = None if not site['details'] else site['details'].split(' ')[-1]
            if ligand_res_num in ligands:
                lig_short = ligands[ligand_res_num]['shortName']
                lig_long = ligands[ligand_res_num]['longName']
                lig_chain = ligands[ligand_res_num]['chainID']
            else:
                ligand_res_num = None
                lig_short = None
                lig_long = None
                lig_chain = None
            for resi in site['site_residues']:
                auth_uniprot_resi_num = resi['author_residue_number']
                ## if residue is the ligand skip it, the keys for %ligands
				## are the 'author residue number' in binding_sites endpoint
                if auth_uniprot_resi_num in ligands:
                    pass
                else:
                    struc_resi_num = resi['residue_number']
                    residue = resi['chem_comp_id']
                    residue_chain = resi['chain_id']
                
                    try:
                        cursor.execute(
                            "REPLACE INTO pdbBindingSites VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
                            (pdb[0], site_ID, struc_resi_num, auth_uniprot_resi_num,
                             residue, residue_chain, ligand_res_num, lig_short,
                             lig_long, lig_chain)
                        )
                        cnx.commit()
                    except mysql.connector.Error as err:
                        print("Something went wrong: {}".format(err))
                        with open(binding_site_log, 'a') as f:
                            f.write(f"{pdb}:   {site_ID}, {auth_uniprot_resi_num}\n")
    n += 1
    print(n)
    #if n == 2000:
        #break
    
