#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 13:18:22 2020

@author: coyote
"""

import neo4j
#print(neo4j.__version__)
import logging
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import numpy as np
import json

"""
Driver for filling Neo4j database using the data queried from MySQL
"""


class ProteinExample:

    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        self.database = None;
    
    
    def close(self):
        self.driver.close()
     
        
    def use_database(self, database):
        self.database=database
    

    # create db function
    def create_db(self, db_name):
        with self.driver.session() as session:
            session.run("CREATE DATABASE $db_name IF NOT EXISTS", db_name=db_name)
            self.database = db_name


    def _drop(self, database):
        with self.driver.session() as session:
            session.run("DROP DATABASE $db_name IF EXISTS", db_name=database)
            

    def create_proteins(self, proteins):
        if not self.database:
            return "Need to call ProteinExample.use_database(database)"
        for protein in proteins:
            self.create_protein(protein)
        

    def create_protein(self, protein):

        with self.driver.session(database=self.database) as session:
    
                result = session.write_transaction(
                    self._create_protein, protein)
                for res in result:
                    if res["merged"]:
                        print("Merged: ", [res["uniprot ID"], 
                                            res["name"],
                                            res["full name"]])
                    else:
                        print("Created: ", [res["uniprot ID"], 
                                            res["name"],
                                            res["full name"]])
    
    
    @staticmethod
    def _create_protein(tx, protein):
        print("PROTEIN:",protein)
        altGeneNames = None if not protein['geneNamesAlternative'] else json.loads(protein['geneNamesAlternative'])
        altProtNames = None if not protein['alternateNames'] else json.loads(protein['alternateNames'])
        
        #Test updating with Merops ID for proteases
        """
        https://neo4j.com/docs/java-reference/current/java-embedded/property-values/
        NULL is not a valid property value. Setting a property to NULL 
        is equivalent to deleting the property. 
        """
        meropsID = protein['meropsID'] if protein['meropsID'] else None
        #meropsID = None if not protein['meropsID'] else protein['meropsID']
        
        # Rare case that geneNamePreferred is blank so use short entry name for Node label
        if protein['geneNamePreferred'] is None:
            protein['geneNamePreferred'] = protein['entryName']
        
        print("*************************************", protein['uniProtID'])
        query = ("""
                 MERGE (a:Protein{name: $name, 
                				  uniprotID: $uniprotID, 
                				  organism: $organism, 
                				  taxid: $taxid})
                    ON CREATE SET a.proteinName = $proteinName,
                                  a.meropsID = $meropsID,
                                  a.altGeneNames = CASE $altGeneNames
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altGeneNames
                                                     END,
                                  a.altProtNames = CASE $altProtNames
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altProtNames
                                                     END,
                                  a.proteinFamily = $proteinFamily,
                                  a.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                	ON MATCH SET  a.meropsID     = CASE a.meropsID
                                                     WHEN null
                                                     THEN $meropsID
                                                     ELSE a.meropsID
                                                     END,
                                  a.proteinName  = CASE a.proteinName
                                                     WHEN null
                                                     THEN $proteinName
                                                     ELSE a.proteinName
                                                     END,                                                    
                                  a.altGeneNames = CASE a.altGeneNames
                                                     WHEN null
                                                     THEN $altGeneNames
                                                     ELSE 
                                                       CASE $altGeneNames
                                                         WHEN null
                                                         THEN a.altGeneNames
                                                         ELSE apoc.coll.toSet(a.altGeneNames + $altGeneNames)
                                                       END
                                                    END,
                                  a.altProtNames = CASE a.altProtNames
                                                     WHEN null
                                                     THEN $altProtNames
                                                     ELSE 
                                                       CASE $altProtNames
                                                         WHEN null
                                                         THEN a.altProtNames
                                                         ELSE apoc.coll.toSet(a.altProtNames + $altProtNames)
                                                       END
                                                    END,
                                  a.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                RETURN a.uniprotID, a.name, a.proteinName, exists(a.lastModified) as onMerge
        """)
        
        result = tx.run(query, name=protein['geneNamePreferred'], 
                        uniprotID=protein['uniProtID'],
                        proteinName=protein['proteinName'], 
                        organism=protein['organism'],
                        taxid=protein['taxid'],
                        proteinFamily=protein['headProteinFamily'],
                        meropsID=protein['meropsID'], 
                        altGeneNames=altGeneNames,
                        altProtNames=altProtNames)
        try:
            return [{"uniprot ID": record["a.uniprotID"], 
                     "name": record["a.name"],
                     "full name": record["a.proteinName"],
                     "merged": record["onMerge"]} 
                    for record in result]
        except ServiceUnavailable as exception:
            logging.error("{} raised an error: \n {}".format(query, exception))
            raise
    
    
    def _delete_protein(self, uniprotID):
        """
        Not used in SubDB+
        """
        if not self.database:
            return "Need to call ProteinExample.use_database"
        with self.driver.session(database=self.database) as session:
            result = session.run("MATCH (n:Protein{uniprotID: $uniprotID}) "
                                 "WITH n "
                                 "DELETE n "
                                 "RETURN n.uniprotID", 
                                 uniprotID=uniprotID)
            
            for res in result:
                return res["n.uniprotID"] 

###############################################################################
# IntAct
###############################################################################
    def create_intact_interactions(self, interactions):
        if not self.database:
            return "Need to call ProteinExample.use_database(database)"
        for interaction in interactions:
            direct = str.split(interaction['direction'], "->")
            ## if self interaction
            if len(np.unique(direct)) == 1:
                with self.driver.session(database=self.database) as session:   
                    result = session.write_transaction(
                        self._create_intact_self_interaction, interaction, direct[0], direct[0])
                for res in result:
                    if res["merged"]:
                        print("Merged: ", [res["from"],
                                           res["interaction"],
                                           res["to"]])
                    else:
                        print("Created: ", [res["from"],
                                            res["interaction"],
                                            res["to"]])
            ## not self interaction       
            else:
                with self.driver.session(database=self.database) as session:
                    result = session.write_transaction(
                            self._create_intact_a_b_interaction, interaction, direct[0], direct[1])
                for res in result:
                    if res["merged"]:
                        print("Merged: ", [res["from"],
                                           res["interaction"],
                                           res["to"]])
                    else:
                        print("Created: ", [res["from"],
                                            res["interaction"],
                                            res["to"]])
                
            
    @staticmethod
    def _create_intact_self_interaction(tx, interaction, from_, to_):
        altGeneNamesA = None if not interaction['altGeneNamesA'] else json.loads(interaction['altGeneNamesA'])
        altProtNamesA = None if not interaction['altProtNamesA'] else json.loads(interaction['altProtNamesA'])
        
        print("*************************************", interaction['interactionID'], interaction['interactorA'])
        
        query = ("""
                 MERGE (a:%s{uniprotID: $interactorA, 
                                   taxid: $taxidA})
                	ON CREATE SET a.name = $nameA,
                				  a.organism = $organismA,
                				  a.altGeneNames = CASE $altGeneNamesA
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altGeneNamesA
                                                     END,
                				  a.altProtNames = CASE $altProtNamesA
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altProtNamesA
                                                     END,
                                  a.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                	ON MATCH SET  a.altGeneNames = CASE a.altGeneNames
                                                     WHEN null
                                                     THEN $altGeneNamesA
                                                     ELSE 
                                                       CASE $altGeneNamesA
                                                         WHEN null
                                                         THEN a.altGeneNames
                                                         ELSE apoc.coll.toSet(a.altGeneNames + $altGeneNamesA)
                                                       END
                                                    END,
                                  a.altProtNames = CASE a.altProtNames
                                                     WHEN null
                                                     THEN $altProtNamesA
                                                     ELSE 
                                                       CASE $altProtNamesA
                                                         WHEN null
                                                         THEN a.altProtNamesA
                                                         ELSE apoc.coll.toSet(a.altProtNames + $altProtNamesA)
                                                       END
                                                    END,
                                  a.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                WITH a, $interactionID AS identifier, $publicationID AS pub, $detectionMethod as method,
                    $intact_null as site, $sourceDB as db
                MERGE (%s)-[i:INTERACTION {name: $interactionType}]->(%s)
                	ON CREATE SET i.entries = [apoc.convert.toSortedJsonMap({
                                                    interactionID:identifier,
                                                    `publicationID(s)`:pub,
                                                    detectionMethod:method,
                                                    `site(s)`:NULL,
                                                    sourceDB:db})],
                                  i.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH  SET i.entries = apoc.coll.toSet(i.entries + 
                                              apoc.convert.toSortedJsonMap({
                                                    interactionID:identifier,
                                                    `publicationID(s)`:pub,
                                                    detectionMethod:method,
                                                    `site(s)`:NULL,
                                                    sourceDB:db})),
                                  i.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')                              
                RETURN %s.uniprotID as from, i.name, %s.uniprotID as to, exists(i.lastModified) as onMerge
        """ % (interaction['nodeTypeA'], from_, to_, from_, to_))
        
        result = tx.run(query, 
                        ## had to use placeholder instead for the node types
                        #nodeTypeA=interaction['nodeTypeA'],
                        interactorA=interaction['interactorA'],
                        taxidA=interaction['taxidA'],
                        nameA=interaction['nameA'],
                        organismA=interaction['organismA'],
                        altGeneNamesA=altGeneNamesA,
                        altProtNamesA=altProtNamesA,
                        interactionID=interaction['interactionID'],
                        publicationID=interaction['publicationID'],
                        detectionMethod=interaction['detectionMethod'],
                        intact_null="IntAct(null)",
                        sourceDB=interaction['sourceDB'],
                        directionFrom=from_, directionTo=to_,
                        interactionType=interaction['interactionType'])
        try:
            return [{"from": record["from"],
                     "interaction": record["i.name"],
                     "to": record["to"],
                     "merged": record["onMerge"]} 
                    for record in result]
        except ServiceUnavailable as exception:
            logging.error("{} raised an error: \n {}".format(query, exception))
            raise
    

    @staticmethod
    def _create_intact_a_b_interaction(tx, interaction, from_, to_):
        altGeneNamesA = None if not interaction['altGeneNamesA'] else json.loads(interaction['altGeneNamesA'])
        altProtNamesA = None if not interaction['altProtNamesA'] else json.loads(interaction['altProtNamesA'])
        altGeneNamesB = None if not interaction['altGeneNamesB'] else json.loads(interaction['altGeneNamesB'])
        altProtNamesB = None if not interaction['altProtNamesB'] else json.loads(interaction['altProtNamesB'])
        
        ## if geneNamePreferred None for A or B
        if not interaction['nameA'] and altGeneNamesA:
            interaction['nameA'] = altGeneNamesA[0]
        elif not interaction['nameA'] and not altGeneNamesA:
            interaction['nameA'] = interaction['interactorA']
            
        if not interaction['nameB'] and altGeneNamesB:
            interaction['nameB'] = altGeneNamesB[0]
        elif not interaction['nameB'] and not altGeneNamesB:
            interaction['nameB'] = interaction['interactorB']
            
        print("*************************************", interaction['interactionID'], interaction['interactorA'])

        query = ("""
                 MERGE (a:%s{uniprotID: $interactorA, 
                                   taxid: $taxidA})
                	ON CREATE SET a.name = $nameA,
                				  a.organism = $organismA,
                				  a.altGeneNames = CASE $altGeneNamesA
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altGeneNamesA
                                                     END,
                				  a.altProtNames = CASE $altProtNamesA
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altProtNamesA
                                                     END,
                                  a.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                	ON MATCH SET  a.altGeneNames = CASE a.altGeneNames
                                                     WHEN null
                                                     THEN $altGeneNamesA
                                                     ELSE 
                                                       CASE $altGeneNamesA
                                                         WHEN null
                                                         THEN a.altGeneNames
                                                         ELSE apoc.coll.toSet(a.altGeneNames + $altGeneNamesA)
                                                       END
                                                    END,
                                  a.altProtNames = CASE a.altProtNamesA
                                                     WHEN null
                                                     THEN $altProtNamesA
                                                     ELSE 
                                                       CASE $altProtNamesA
                                                         WHEN null
                                                         THEN a.altProtNamesA
                                                         ELSE apoc.coll.toSet(a.altProtNamesA + $altProtNamesA)
                                                       END
                                                    END,
                                  a.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')

                 MERGE (b:%s{uniprotID: $interactorB, 
                                   taxid: $taxidB})
                	ON CREATE SET b.name = $nameB,
                				  b.organism = $organismB,
                				  b.altGeneNames = CASE $altGeneNamesB
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altGeneNamesB
                                                     END,
                				  b.altProtNames = CASE $altProtNamesB
                                                     WHEN null
                                                     THEN NULL
                                                     ELSE $altProtNamesB
                                                     END,
                                  b.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                	ON MATCH SET  b.altGeneNames = CASE b.altGeneNames
                                                     WHEN null
                                                     THEN $altGeneNamesB
                                                     ELSE 
                                                       CASE $altGeneNamesB
                                                         WHEN null
                                                         THEN b.altGeneNames
                                                         ELSE apoc.coll.toSet(b.altGeneNames + $altGeneNamesB)
                                                       END
                                                    END,
                                  b.altProtNames = CASE b.altProtNamesB
                                                     WHEN null
                                                     THEN $altProtNamesB
                                                     ELSE 
                                                       CASE $altProtNamesB
                                                         WHEN null
                                                         THEN b.altProtNamesB
                                                         ELSE apoc.coll.toSet(b.altProtNamesB + $altProtNamesB)
                                                       END
                                                    END,
                                  b.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                 
                
                WITH a, b, $interactionID AS identifier, $publicationID AS pub, $detectionMethod as method,
                    $intact_null as site, $sourceDB as db
                MERGE (%s)-[i:INTERACTION {name: $interactionType}]->(%s)
                	ON CREATE SET i.entries = [apoc.convert.toSortedJsonMap({
                                                    interactionID:identifier,
                                                    `publicationID(s)`:pub,
                                                    detectionMethod:method,
                                                    `site(s)`:NULL,
                                                    sourceDB:db})],
                                  i.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH  SET i.entries = apoc.coll.toSet(i.entries + 
                                              apoc.convert.toSortedJsonMap({
                                                    interactionID:identifier,
                                                    `publicationID(s)`:pub,
                                                    detectionMethod:method,
                                                    `site(s)`:NULL,
                                                    sourceDB:db})),
                                  i.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')                              
                RETURN %s.uniprotID as from, i.name, %s.uniprotID as to, exists(i.lastModified) as onMerge
        """ % (interaction['nodeTypeA'], interaction['nodeTypeB'], from_, to_, from_, to_))
        
        result = tx.run(query, 
                        interactorA=interaction['interactorA'],
                        taxidA=interaction['taxidA'],
                        nameA=interaction['nameA'],
                        organismA=interaction['organismA'],
                        altGeneNamesA=altGeneNamesA,
                        altProtNamesA=altProtNamesA,
                        interactorB=interaction['interactorB'],
                        taxidB=interaction['taxidB'],
                        nameB=interaction['nameB'],
                        organismB=interaction['organismB'],
                        altGeneNamesB=altGeneNamesB,
                        altProtNamesB=altProtNamesB,
                        interactionID=interaction['interactionID'],
                        publicationID=interaction['publicationID'],
                        detectionMethod=interaction['detectionMethod'],
                        intact_null="IntAct(null)",
                        sourceDB=interaction['sourceDB'],
                        directionFrom=from_, directionTo=to_,
                        interactionType=interaction['interactionType'])
        try:
            return [{"from": record["from"],
                     "interaction": record["i.name"],
                     "to": record["to"],
                     "merged": record["onMerge"]} 
                    for record in result]
        except ServiceUnavailable as exception:
            logging.error("{} raised an error: \n {}".format(query, exception))
            raise    
###############################################################################    
# MEROPS
###############################################################################
    def create_merops_interactions(self, interactions):
        if not self.database:
            return "Need to call ProteinExample.use_database(database)"
        for interaction in interactions:
            with self.driver.session(database=self.database) as session:   
                 result = session.write_transaction(
                     self._create_merops_a_b_interaction, interaction)
            for res in result:
                if res["merged"]:
                    print("Merged: ", [res["from"],
                                       res["interaction"],
                                       res["to"]])
                else:
                    print("Created: ", [res["from"],
                                        res["interaction"],
                                        res["to"]])
                   

    @staticmethod
    def _create_merops_a_b_interaction(tx, interaction):
        subGeneAlt = None if not interaction['subGeneAlt'] else json.loads(interaction['subGeneAlt'])

        # Only doing proteases from UniProt that have merops IDs
        # so will always be in data warehouse initially
        # if we cannot get subsrate info from Uniprot make the 
        # Uniprot id from Merops database the node name
        query = ("""
                 MERGE (a:Protein{uniprotID: $proteaseUniprot, 
                                   taxid: $proteaseTaxId})
                    ON CREATE SET a.name = $proteaseGenePreferred,
                                  a.meropsID = $meropsID,
                                  a.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH SET  a.meropsID = $meropsID,
                                  a.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                 
                 MERGE (b:Protein{uniprotID: $substrateUniprot})
                	ON CREATE SET b.name = CASE $substrateGenePreferred
                                             WHEN null
                                             THEN $substrateUniprot
                                             ELSE $substrateGenePreferred
                                             END,
                                  b.proteinName = $substrateName, 
                                  b.altGeneNames = $subGeneAlt,
                				  b.organism = $substrateOrganism,
                                  b.taxid = CASE $substrateTax
                                              WHEN null
                                              THEN $proteaseTaxId
                                              ELSE $substrateTax
                                              END,
                                  b.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH SET  b.proteinName  = CASE b.proteinName
                                                     WHEN null
                                                     THEN $substrateName
                                                     ELSE b.proteinName
                                                     END,  
                                  b.altGeneNames = CASE b.altGeneNames
                                                     WHEN null
                                                     THEN $subGeneAlt
                                                     ELSE 
                                                         CASE $subGeneAlt
                                                           WHEN null
                                                           THEN b.altGeneNames
                                                         ELSE apoc.coll.toSet(b.altGeneNames + $subGeneAlt)
                                                         END
                                                   END,
                                  b.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                 
                WITH a, b, $cleavageID AS identifier, $Ref AS pub, $cleavage_type as method,
                    $Substrate_formula as subsite, 'MEROPS' as db
                MERGE (a)-[i:INTERACTION {name: $interactionType}]->(b)
                	ON CREATE SET i.entries = [apoc.convert.toSortedJsonMap({
                                                    interactionID:identifier,
                                                    `publicationID(s)`:pub,
                                                    detectionMethod:method,
                                                    `site(s)`:subsite,
                                                    sourceDB:db})],
                                  i.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH  SET i.entries = apoc.coll.toSet(i.entries + 
                                              apoc.convert.toSortedJsonMap({
                                                    interactionID:identifier,
                                                    `publicationID(s)`:pub,
                                                    detectionMethod:method,
                                                    `site(s)`:subsite,
                                                    sourceDB:db})),
                                  i.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')                              
                RETURN a.uniprotID as from, i.name, b.uniprotID as to, exists(i.lastModified) as onMerge
        """)
        
        result = tx.run(query, 
                        proteaseUniprot=interaction['proteaseUniprot'],
                        proteaseTaxId=interaction['proteaseTaxId'],
                        proteaseGenePreferred=interaction['proteaseGenePreferred'],
                        meropsID=interaction['meropsID'],
                        substrateUniprot=interaction['substrateUniprot'],
                        substrateGenePreferred=interaction['substrateGenePreferred'],
                        substrateName=interaction['substrateName'],
                        subGeneAlt=subGeneAlt,
                        substrateOrganism=interaction['substrateOrganism'],
                        substrateTax=interaction['substrateTax'],
                        cleavageID=interaction['cleavageID'],
                        Ref=interaction['Ref'],
                        cleavage_type=interaction['cleavage_type'],
                        Substrate_formula=interaction['Substrate_formula'],
                        interactionType="proteolysis")
        try:
            return [{"from": record["from"],
                     "interaction": record["i.name"],
                     "to": record["to"],
                     "merged": record["onMerge"]} 
                    for record in result]
        except ServiceUnavailable as exception:
            logging.error("{} raised an error: \n {}".format(query, exception))
            raise       
    
    
###############################################################################
# PSP
###############################################################################    
    def create_psp_interactions(self, interactions):
        if not self.database:
            return "Need to call ProteinExample.use_database(database)"
        for interaction in interactions:
            with self.driver.session(database=self.database) as session:   
                 result = session.write_transaction(
                     self._create_psp_a_b_interaction, interaction)
            for res in result:
                if res["merged"]:
                    print("Merged: ", [res["from"],
                                       res["interaction"],
                                       res["to"]])
                else:
                    print("Created: ", [res["from"],
                                        res["interaction"],
                                        res["to"]])
                   

    @staticmethod
    def _create_psp_a_b_interaction(tx, interaction):
        geneNameAltSub = None if not interaction['geneNameAltSub'] else [interaction['geneNameAltSub']]

        if interaction['inVivo'] and interaction['inVitro']:
            method = "in Vivo, in Vitro"
        elif interaction['inVivo'] and not interaction['inVitro']:
            method = "in Vivo"
        elif interaction['inVitro'] and not interaction['inVivo']:
            method = "in Vitro"    
        else:
            method = None
            
            
        query = ("""
                 MERGE (a:Protein{uniprotID: $uniProtIDKin, 
                                   taxid: $kinTaxid})
                    ON CREATE SET a.name = $geneNamePreferredKin,
                                  a.proteinName = $protNamePreferredKin,
                                  a.organism = $kinOrganism,
                                  a.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH SET  a.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                 
                 MERGE (b:Protein{uniprotID: $uniProtIDSub,
                                   taxid: $subTaxid})
                	ON CREATE SET b.name = $geneNamePreferredSub,
                                  b.proteinName = $protNamePreferredSub,
                                  b.organism = $subOrganism,
                                  b.altGeneNames = $geneNameAltSub,            
                                  b.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH SET  b.altGeneNames = CASE b.altGeneNames
                                                     WHEN null
                                                     THEN $geneNameAltSub
                                                     ELSE 
                                                       CASE $geneNameAltSub
                                                         WHEN null
                                                         THEN b.altGeneNames
                                                         ELSE apoc.coll.toSet(b.altGeneNames + $geneNameAltSub)
                                                       END
                                                    END,
                                  b.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                 
                WITH a, b, $sitePlusMinus7AA AS identifier, $psp_null AS pub, $method as method,
                    $subModSite as subsite, 'PhosphoSitePlus' as db
                MERGE (a)-[i:INTERACTION {name: $interactionType}]->(b)
                	ON CREATE SET i.entries = [apoc.convert.toSortedJsonMap({
                                                    sitePlusMinus7AA:identifier,
                                                    `publicationID(s)`:NULL,
                                                    detectionMethod:method,
                                                    `site(s)`:subsite,
                                                    sourceDB:db})],
                                  i.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH  SET i.entries = apoc.coll.toSet(i.entries + 
                                              apoc.convert.toSortedJsonMap({
                                                    sitePlusMinus7AA:identifier,
                                                    `publicationID(s)`:NULL,
                                                    detectionMethod:method,
                                                    `site(s)`:subsite,
                                                    sourceDB:db})),
                                  i.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')                              
                RETURN a.uniprotID as from, i.name, b.uniprotID as to, exists(i.lastModified) as onMerge
        """)
        
        result = tx.run(query, 
                        uniProtIDKin=interaction['uniProtIDKin'],
                        kinTaxid=interaction['kinTaxid'],
                        geneNamePreferredKin=interaction['geneNamePreferredKin'],
                        protNamePreferredKin=interaction['protNamePreferredKin'],
                        kinOrganism=interaction['kinOrganism'],
                        uniProtIDSub=interaction['uniProtIDSub'],
                        subTaxid=interaction['subTaxid'],
                        geneNamePreferredSub=interaction['geneNamePreferredSub'],
                        protNamePreferredSub=interaction['protNamePreferredSub'],
                        subOrganism=interaction['subOrganism'],
                        geneNameAltSub=geneNameAltSub,
                        sitePlusMinus7AA=interaction['sitePlusMinus7AA'],
                        psp_null="PSP(null)",
                        method=method,
                        subModSite=interaction['subModSite'],
                        interactionType="phosphorylation")
        try:
            return [{"from": record["from"],
                     "interaction": record["i.name"],
                     "to": record["to"],
                     "merged": record["onMerge"]} 
                    for record in result]
        except ServiceUnavailable as exception:
            logging.error("{} raised an error: \n {}".format(query, exception))
            raise



###############################################################################
# DEPOD
###############################################################################
    def create_depod_interactions(self, interactions):
        if not self.database:
            return "Need to call ProteinExample.use_database(database)"
        for interaction in interactions:
            with self.driver.session(database=self.database) as session:   
                 result = session.write_transaction(
                     self._create_depod_a_b_interaction, interaction)
            for res in result:
                if res["merged"]:
                    print("Merged: ", [res["from"],
                                       res["interaction"],
                                       res["to"]])
                else:
                    print("Created: ", [res["from"],
                                        res["interaction"],
                                        res["to"]])
                   

    @staticmethod
    def _create_depod_a_b_interaction(tx, interaction):

        if interaction['inVivo'] and interaction['inVitro']:
            method = "in Vivo, in Vitro"
        elif interaction['inVivo'] and not ['inVitro']:
            method = "in Vivo"
        elif interaction['inVitro'] and not ['inVivo']:
            method = "in Vitro"    
        else:
            method = None
            
            
        query = ("""
                 MERGE (a:Protein{uniprotID: $uniProtIDPPase, 
                                   taxid: $ppaseTaxid})
                    ON CREATE SET a.name = $geneNamePreferredPPase,
                                  a.proteinName = $protNamePreferredPPase,
                                  a.organism = $ppaseOrganism,
                                  a.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH SET  a.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                 
                 MERGE (b:Protein{uniprotID: $uniProtIDSub,
                                   taxid: $subTaxid})
                	ON CREATE SET b.name = $geneNamePreferredSub,
                                  b.proteinName = $protNamePreferredSub,
                                  b.organism = $subOrganism,            
                                  b.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH SET  b.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                 
                WITH a, b, $sitePlusMinus5AA AS identifier, $literature AS pub, $method as method,
                    $subDephospoSites as subsite, 'DEPhOsphorylation Database' as db
                MERGE (a)-[i:INTERACTION {name: $interactionType}]->(b)
                	ON CREATE SET i.entries = [apoc.convert.toSortedJsonMap({
                                                    sitePlusMinus5AA:identifier,
                                                    `publicationID(s)`:pub,
                                                    detectionMethod:method,
                                                    `site(s)`:subsite,
                                                    sourceDB:db})],
                                  i.created = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')
                    ON MATCH  SET i.entries = apoc.coll.toSet(i.entries + 
                                              apoc.convert.toSortedJsonMap({
                                                    sitePlusMinus5AA:identifier,
                                                    `publicationID(s)`:NULL,
                                                    detectionMethod:method,
                                                    `site(s)`:subsite,
                                                    sourceDB:db})),
                                  i.lastModified = apoc.date.format(timestamp(),'ms','yyyy-MM-dd HH:mm:ss.sss','EST')                              
                RETURN a.uniprotID as from, i.name, b.uniprotID as to, exists(i.lastModified) as onMerge
        """)
        
        result = tx.run(query, 
                        uniProtIDPPase=interaction['uniProtIDPPase'],
                        ppaseTaxid=interaction['ppaseTaxid'],
                        geneNamePreferredPPase=interaction['geneNamePreferredPPase'],
                        protNamePreferredPPase=interaction['protNamePreferredPPase'],
                        ppaseOrganism=interaction['ppaseOrganism'],
                        uniProtIDSub=interaction['uniProtIDSub'],
                        subTaxid=interaction['subTaxid'],
                        geneNamePreferredSub=interaction['geneNamePreferredSub'],
                        protNamePreferredSub=interaction['protNamePreferredSub'],
                        subOrganism=interaction['subOrganism'],
                        sitePlusMinus5AA=interaction['sitePlusMinus5AA'],
                        literature=interaction['literature'],
                        method=method,
                        subDephospoSites=interaction['subDephospoSites'],
                        interactionType="dephosphorylation")
        try:
            return [{"from": record["from"],
                     "interaction": record["i.name"],
                     "to": record["to"],
                     "merged": record["onMerge"]} 
                    for record in result]
        except ServiceUnavailable as exception:
            logging.error("{} raised an error: \n {}".format(query, exception))
            raise



###############################################################################
    # TODO: delete interactions
    # Not Implemented
    def _delete_interaction(self, db, prot_a, up_or_down, prot_b):
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
    
    # Not Implemented
    def _get_recursive_n_interactions(self, protein, n_edges):
        """
        Notes
        -------
        Not Implemented \n
        Implemented in Shiny App
        """   
        with self.driver.session(database=self.database) as session:
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
          

