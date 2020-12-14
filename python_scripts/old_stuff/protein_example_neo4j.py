#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 13:18:22 2020

@author: coyote
"""


import logging
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable


class ProteinExample:

    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
    
    
    def close(self):
        self.driver.close()
    
    """
    Get all nodes
    MATCH (n)-[r]->(m)
    RETURN n,r,m
    """

    # create db function
    def create_db(self, db_name):
        with self.driver.session() as session:
            session.run("CREATE DATABASE $db_name IF NOT EXISTS", db_name=db_name)


    def drop(self, db_name):
        with self.driver.session() as session:
            session.run("DROP DATABASE $db_name IF EXISTS", db_name=db_name)
            

    def create_proteins(self, database, prot_name_type_pairs):
        """
        
        Parameters
        ----------
        database : TYPE
            DESCRIPTION.
        prot_name_type_pairs : Collection
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for prot, type_ in prot_name_type_pairs:
            self.create_protein(database, prot, type_)
        

    def create_protein(self, database, prot_name, type_):
        """
        
        Parameters
        ----------
        database : TYPE
            DESCRIPTION.
        prot_name : TYPE
            DESCRIPTION.
        type_ : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """  

        with self.driver.session(database=database) as session:
    
                result = session.write_transaction(
                    self._create_protein, prot_name=prot_name, type_=type_)
                for res in result:
                    print("Created: ", [res["protein"], res["type"]])
    
    
    @staticmethod
    def _create_protein(tx, prot_name, type_):
        query = (
            "MERGE (n:Protein { name: $prot_name, classification: $type_ }) "
            "RETURN n.name, n.classification"
        )
        result = tx.run(query, prot_name=prot_name, type_=type_)
        try:
            return [{"protein": record["n.name"], "type": record["n.classification"]} 
                    for record in result]
        except ServiceUnavailable as exception:
            logging.error("{} raised an error: \n {}".format(query, exception))
            raise

    
    def create_interactions(self, database, interaction_set):
        """
        
        Parameters
        ----------
        database : TYPE
            DESCRIPTION.
        interaction_set : Collection
            DESCRIPTION.

        Returns
        -------
        None.

        """
    
        for a_name, direct, b_name in interaction_set:
            self.create_interaction(database, a_name, direct, b_name)


    def create_interaction(self, database, prot_a, up_or_down, prot_b):
        if up_or_down not in ["+", "-"]:
            return "up_or_down only takes '+' or '-'"
        
        with self.driver.session(database=database) as session:
            result = session.run("MATCH (a:Protein{name: $prot_a}) "
                                 "MATCH (b:Protein{name: $prot_b}) "
                                 "MERGE (a)-[r:REGULATES {direction:$up_or_down}]->(b) "
                                 "RETURN a.name, r.direction, b.name", 
                                 prot_a=prot_a, up_or_down=up_or_down, prot_b=prot_b)  
            for res in result:
                print([res["a.name"], res["r.direction"], res["b.name"]])

    
    def delete_interaction(self, db, prot_a, up_or_down, prot_b):
        if up_or_down not in ["+", "-"]:
            return "up_or_down only takes '+' or '-'"
        
        with self.driver.session(database=db) as session:
            result = session.run("MATCH (a:Protein),(b:Protein) "
                                 "WHERE a.name = $prot_a AND b.name = $prot_b "
                                 "MATCH (a)-[r:REGULATES {direction:$up_or_down}]->(b) "
                                 "WITH r, a, b, r.direction AS direction "
                                 "DELETE r "
                                 "RETURN a.name, direction, b.name", 
                                 prot_a=prot_a, up_or_down=up_or_down, prot_b=prot_b)  
            
            for res in result:
                return [res["a.name"], res["direction"], res["b.name"]] 
            
          
    def get_recursive_n_interactions(self, db, protein, n_edges):
        """

        Parameters
        ----------
        db : TYPE
            DESCRIPTION.
        protein : TYPE
            DESCRIPTION.
        n_edges : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

     
        MATCH p = (n:Protein { name:'C' })<-[:REGULATES*1..2]->(b:Protein) 
        WITH *, relationships(p) as rs
        RETURN startNode(last(rs)).name as Protein1, 
            last(rs).direction as Regulates, 
            endNode(last(rs)).name as Protein2, 
            length(p)
        """   
        with self.driver.session(database=db) as session:
            result = session.run("MATCH p = (n:Protein { name:$prot })<-[:REGULATES*1..%d]->(b:Protein) \
                                  WITH *, relationships(p) AS rs \
                                  RETURN \
                                     startNode(last(rs)).name AS Protein1, \
                                     last(rs).direction AS Regulates, \
                                     endNode(last(rs)).name AS Protein2, \
                                     length(p) AS pathLength" % n_edges,
                                {"prot": protein})
            for res in result:
                print([res["Protein1"], res["Regulates"], res["Protein2"], res["pathLength"]])
          



if __name__ == "__main__":
    
    # desktop version
    #uri = "neo4j://localhost:11003"
    # some times needs this
    uri = "neo4j://localhost:7687"
    user = "neo4j"
    # you'll need a password after you start 
    password = "********"
    database = "protein-test"
  
    prot_db = ProteinExample(uri, user, password)
        
    prot_db.drop(database)
    prot_db.create_db(database)
    
    proteins = [('A', 'kinase'),
                ('B', 'phosphatase'),
                ('C', 'kinase'),
                ('D', 'kinase'),
                ('E', 'kinase'),
                ('F', 'phosphatase'),
                ('G', 'transferase'),
                ('H', 'transferase'),
                ('T', 'methyl_transferase'),
                ('M', 'methyl_transferase'),
                ('L', 'phosphatase')]
    prot_db.create_proteins(database, proteins)
    
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
    prot_db.create_interactions(database, regulations)
    
    # delete
    del_regulation = ('A', '-', 'B')
    prot_db.delete_interaction(database, del_regulation[0], del_regulation[1], del_regulation[2])
    del_regulation = ('L', '+', 'A')
    prot_db.delete_interaction(database, del_regulation[0], del_regulation[1], del_regulation[2])
    
    # re add
    regulation =  ('A', '-', 'B')
    prot_db.create_interaction(database, regulation[0], regulation[1], regulation[2])     
    regulation =  ('L', '+', 'A')
    prot_db.create_interaction(database, regulation[0], regulation[1], regulation[2])  
    
    # The whole kit and kabuttle
    print("\nHere is what you have been waiting for\n")
    prot_db.get_recursive_n_interactions(database, 'C', 2)
