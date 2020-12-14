#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 03:16:37 2020

@author: coyote
"""

from neo4j import GraphDatabase
import neo4j
print(neo4j.__version__)

import MySQLdb

# community version
#uri = "neo4j://localhost:7687"
#driver = GraphDatabase.driver(uri, auth=("neo4j", "Swimgolf1212**"))

# desktop version
uri = "neo4j://localhost:11003"
uri = "bolt://localhost:7687"
driver = GraphDatabase.driver(uri, auth=("neo4j", "Swimgolf1212**"))
#driver = GraphDatabase.driver(uri, auth=("neo4j", "password"))


"""
Get all nodes
MATCH (n)-[r]->(m)
RETURN n,r,m
"""

# create db function
def create_db(driver, db_name):
    with driver.session() as session:
        session.run("CREATE DATABASE $db_name IF NOT EXISTS", db_name=db_name)
    driver.close() 

# create db
create_db(driver, "protein-test")

    
# add nodes function
def create_protein(tx, prot_name, type_):
    result = tx.run("CREATE (n:Protein { name: $prot_name, classification: $type_ }) "
                        "RETURN n.name, n.classification",
                    prot_name=prot_name, type_=type_)
    
    for res in result:
        return [res["n.name"], res["n.classification"]]
    
# test
with driver.session(database="protein-test") as session:
    prot = session.write_transaction(create_protein, 'A', 'kinase')
    print(prot) 
driver.close()  

proteins = [('B', 'phosphatase'),
            ('C', 'peptidase'),
            ('D', 'isomerase'),
            ('E', 'kinase'),
            ('F', 'phoshorylase'),
            ('G', 'ubiquitinase'),
            ('H', 'transferase'),
            ('T', 'methyl transferase'),
            ('M', 'acetyl transferase'),
            ('L', 'kinase')]

with driver.session(database="protein-test") as session:
    for prot in proteins:
        entry = session.write_transaction(create_protein, prot[0], prot[1])
        print(entry)         
driver.close()  


# interaction function
def create_interaction(tx, prot_a, up_or_down, prot_b):
    if up_or_down not in ["+", "-"]:
        return "up_or_down only takes '+' or '-'"
    result = tx.run("MATCH (a:Protein),(b:Protein) "
                    "WHERE a.name = $prot_a AND b.name = $prot_b "
                    "MERGE (a)-[r:REGULATES {direction:$up_or_down}]->(b) "
                    "RETURN a.name, r.direction, b.name", 
                    prot_a=prot_a, up_or_down=up_or_down, prot_b=prot_b)  
    
    for res in result:
        return [res["a.name"], res["r.direction"], res["b.name"]] 

regulations = [('A', '-', 'B'),
               ('B', '-', 'C'),
               ('C', '+', 'D'),
               ('D', '-', 'E'),
               ('E', '+', 'F'),
               ('F', '-', 'G'),
               ('H', '-', 'T'),
               ('T', '-', 'M'),
               ('M', '+', 'C'),
               ('C', '-', 'L'),
               ('L', '+', 'A')]

with driver.session(database="protein-test") as session:
    for reg in regulations:
        regulation = session.write_transaction(create_interaction, 
                                               reg[0], reg[1], reg[2])
        print(regulation)       
driver.close()  


# remove interaction function
def delete_interaction(tx, prot_a, up_or_down, prot_b):
    if up_or_down not in ["+", "-"]:
        return "up_or_down only takes '+' or '-'"
    result = tx.run("MATCH (a:Protein),(b:Protein) "
                    "WHERE a.name = $prot_a AND b.name = $prot_b "
                    "MATCH (a)-[r:REGULATES {direction:$up_or_down}]->(b) "
                    "WITH r, a, b, r.direction AS direction "
                    "DELETE r "
                    "RETURN a.name, direction, b.name", 
                    prot_a=prot_a, up_or_down=up_or_down, prot_b=prot_b)  
    
    for res in result:
        return [res["a.name"], res["direction"], res["b.name"]] 

del_regulation = ('A', '-', 'B')
readd_regulation = ('A', '-', 'B')

# test delete
with driver.session(database="protein-test") as session:
    deleted = session.write_transaction(delete_interaction, 
                                        del_regulation[0], 
                                        del_regulation[1], 
                                        del_regulation[2])
    print(deleted)        
driver.close()  


# add back
with driver.session(database="protein-test") as session:
    readded = session.write_transaction(create_interaction, 
                                        readd_regulation[0], 
                                        readd_regulation[1], 
                                        readd_regulation[2])
    print(readded) 
driver.close()  
        


"""
MATCH p = (n:Protein { name:'C' })<-[:REGULATES*1..2]->(b:Protein) 
WITH *, relationships(p) as rs
RETURN startNode(last(rs)).name as Protein1, 
    last(rs).direction as Regulates, 
    endNode(last(rs)).name as Protein2, 
    length(p)
"""       
#@staticmethod        
def get_recursive_n_interactions(tx, protein, n_edges):
    result = tx.run("MATCH p = (n:Protein { name:$prot })<-[:REGULATES*1..%d]->(b:Protein) \
                     WITH *, relationships(p) AS rs \
                     RETURN \
                         startNode(last(rs)).name AS Protein1, \
                         last(rs).direction AS Regulates, \
                         endNode(last(rs)).name AS Protein2, \
                         length(p) AS pathLength" % n_edges,
                    {"prot": protein})
    for res in result:
        print([res["Protein1"], res["Regulates"], res["Protein2"], res["pathLength"]])
      
with driver.session(database="protein-test") as session:
    session.read_transaction(get_recursive_n_interactions, 'C', 2)
driver.close()
    
