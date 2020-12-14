#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 02:11:04 2020

@author: coyote
"""


##############################################################################
# Test queries for neo4j
##############################################################################

# testing getting entries that overlap in merops and intact to test neo4j
query_merops = """
(SELECT  DISTINCT m.cleavageID, u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, m.proteaseUniprot
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorA AND i.interactorB IN (SELECT u.uniProtID FROM proteinsUniprot u)
JOIN 
	( SELECT sel.cleavageID, sel.Substrate_name as substrateName, prot.geneNamePreferred as substrateGenePreferred,
			sel.substrateUniprot, sel.substrateOrganism, 
			sel.protease, sel.`code`, sel.proteaseUniprot, sel.geneNamePreferred as proteaseGenePreferred,
			sel.proteinName as proteaseName, sel.alternateNames as proteaseAltNames, 
			sel.taxid as proteaseTaxId, sel.organism as proteaseOrganism
		FROM
		(SELECT s.cleavageID, s.Substrate_name, s.Uniprot AS substrateUniprot, s.organism as substrateOrganism, 
			s.protease, s.`code`, p.uniProtID as proteaseUniprot, p.geneNamePreferred, p.proteinName, 
			p.alternateNames, p.taxid, p.organism
		FROM protTest.Substrate_search s
		INNER JOIN 
			proteinsUniprot p
			ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%', p.organism, '%')
		WHERE p.uniProtID IN
			( SELECT i.interactorA
			  FROM intact i 
			  GROUP BY i.interactorA)
		) as sel
		LEFT JOIN
			proteinsUniprot prot
			ON sel.substrateUniprot = prot.uniProtID 
	) as m
    ON u.uniProtID = m.proteaseUniprot
LIMIT 10)

UNION ALL

(SELECT DISTINCT m.cleavageID, u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, m.proteaseUniprot
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorB AND i.interactorA IN (SELECT u.uniProtID FROM proteinsUniprot u)
JOIN 
	( SELECT sel.cleavageID, sel.Substrate_name as substrateName, prot.geneNamePreferred as substrateGenePreferred,
			sel.substrateUniprot, sel.substrateOrganism, 
			sel.protease, sel.`code`, sel.proteaseUniprot, sel.geneNamePreferred as proteaseGenePreferred,
			sel.proteinName as proteaseName, sel.alternateNames as proteaseAltNames, 
			sel.taxid as proteaseTaxId, sel.organism as proteaseOrganism
		FROM
		(SELECT s.cleavageID, s.Substrate_name, s.Uniprot AS substrateUniprot, s.organism as substrateOrganism, 
			s.protease, s.`code`, p.uniProtID as proteaseUniprot, p.geneNamePreferred, p.proteinName, 
			p.alternateNames, p.taxid, p.organism
		FROM protTest.Substrate_search s
		INNER JOIN 
			proteinsUniprot p
			ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%', p.organism, '%')
		WHERE p.uniProtID IN
			( SELECT i.interactorB
			  FROM intact i 
			  GROUP BY i.interactorB)
		) as sel
		LEFT JOIN
			proteinsUniprot prot
			ON sel.substrateUniprot = prot.uniProtID 
	) as m
   ON u.uniProtID = m.proteaseUniprot
LIMIT 10);
"""

# testing getting entries that overlap in phospositeplus and intact to test neo4j
query_psp = """
(SELECT DISTINCT u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, k.uniProtIDKin
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorA AND i.interactorB IN (SELECT u.uniProtID FROM proteinsUniprot u)
JOIN 
	( SELECT DISTINCT uniProtIDKin
      FROM kinasePhosphoSitePlus 
	) as k
    ON u.uniProtID = k.uniProtIDKin
LIMIT 10)

UNION ALL

(SELECT DISTINCT u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, k.uniProtIDKin
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorB AND i.interactorA IN (SELECT u.uniProtID FROM proteinsUniprot u)
JOIN 
	( SELECT DISTINCT uniProtIDKin
      FROM kinasePhosphoSitePlus 
	) as k
    ON u.uniProtID = k.uniProtIDKin
LIMIT 10);
"""

# testing getting entries that overlap in DEPOD and intact to test neo4j
query_ppase = """
(SELECT DISTINCT u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, d.uniProtIDPPase
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorA AND i.interactorB IN (SELECT u.uniProtID FROM proteinsUniprot u)
JOIN 
	( SELECT DISTINCT d.uniProtIDPPase
      FROM depodPhosphatase d 
	) as d
    ON u.uniProtID = d.uniProtIDPPase
LIMIT 10) 

UNION ALL

(SELECT DISTINCT u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, d.uniProtIDPPase
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorB AND i.interactorA IN (SELECT u.uniProtID FROM proteinsUniprot u)
JOIN 
	( SELECT DISTINCT d.uniProtIDPPase
      FROM depodPhosphatase d 
	) as d
    ON u.uniProtID = d.uniProtIDPPase
LIMIT 10);
"""


##############################################################################
# Final queries for neo4j
# Commenting out lines used for testing
##############################################################################

# query uniprot to initially load neo4j
query_up = """
    SELECT u.geneNamePreferred, u.uniProtID, u.proteinName, u.organism, u.entryName,
    	u.taxid, u.geneNamesAlternative, u.alternateNames, u.headProteinFamily, u.meropsID
    FROM proteinsUniprot u
    /*WHERE u.geneNamesAlternative IS NULL*/
    /*WHERE u.uniProtID IN (%s)*/;
"""

# query intact to initially load neo4j
query_intact = """
SELECT i.interactionID, i.interactorA, i.interactorB, i.taxidA, i.taxidB,
	i.geneNamePreferredA as nameA, i.geneNamePreferredB as nameB, i.organismA, i.organismB,
    i.alternateNamesA as altProtNamesA, i.alternateNamesB as altProtNamesB, 
    i.geneNamesAlternativeA as altGeneNamesA, i.geneNamesAlternativeB as altGeneNamesB,
    i.publicationID, i.detectionMethod, i.interactionType, sourceDB,
    CASE
		WHEN i.interactorA REGEXP "^EBI-" THEN "Molecule"
        WHEN i.interactorA REGEXP "^[0-9]" THEN "Molecule"
        ELSE "Protein"
	END as nodeTypeA,
    CASE
		WHEN i.interactorB REGEXP "^EBI-" THEN "Molecule"
        WHEN i.interactorB REGEXP "^[0-9]" THEN "Molecule"
        ELSE "Protein"
	END as nodeTypeB, lower(i.direction) as direction
FROM protTest.intact i
/*WHERE interactionID IN (%s)*/
"""

# query merops to initially load neo4j
# SEE merops_view.sql in SQL folder
query_merops2 = """
SELECT sel.cleavageID, sel.Ref, sel.cleavage_type, sel.Substrate_formula,
	sel.Substrate_name as substrateName, m.geneNamePreferred as substrateGenePreferred,
    m.geneNamesAlternative as subGeneAlt, m.taxid as substrateTax,
	sel.substrateUniprot, sel.substrateOrganism, 
	sel.protease, sel.`code`, sel.proteaseUniprot, sel.geneNamePreferred as proteaseGenePreferred,
    sel.proteinName as proteaseName, sel.alternateNames as proteaseAltNames, 
    sel.taxid as proteaseTaxId, sel.organism as proteaseOrganism, sel.meropsID
FROM
(SELECT s.cleavageID, s.Ref, s.cleavage_type, s.Substrate_formula,
	s.Substrate_name, SUBSTRING(s.Uniprot, 1, 6) AS substrateUniprot, s.organism as substrateOrganism, 
	s.protease, s.`code`, p.uniProtID as proteaseUniprot, p.geneNamePreferred, p.proteinName, 
    p.alternateNames, p.taxid, p.organism, p.meropsID
FROM protTest.Substrate_search s
INNER JOIN 
	protTest.proteinsUniprot p
    ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%%', p.organism, '%%')) as sel
LEFT JOIN
	meropsSubsUniprot m
    ON sel.substrateUniprot = m.uniProtID
/*WHERE sel.cleavageID IN (%s) AND sel.substrateUniprot IS NOT NULL*/
WHERE sel.substrateUniprot IS NOT NULL
"""

# query phosphositeplus to initially load neo4j
query_psp2 = """
SELECT uniProtIDKin, geneNamePreferredKin, protNamePreferredKin, kinTaxid, kinOrganism, 
	uniProtIDSub, geneNamePreferredSub, geneNameAltSub, protNamePreferredSub, subTaxid, subOrganism,
    subModSite, inVivo, inVitro, sitePlusMinus7AA
FROM kinasePhosphoSitePlus
/*WHERE uniProtIDSub IN (%s)
LIMIT 15*/;
"""

# query DEPOD ppases to initially load neo4j
query_ppase2 = """
SELECT d.uniProtIDPPase, d.geneNamePreferredPPase, protNamePreferredPPase, d.ppaseTaxid, d.ppaseOrganism,
	d.uniProtIDSub, d.geneNamePreferredSub, protNamePreferredSub, d.subTaxid, d.subOrganism,
    d.sitePlusMinus5AA, d.literature, d.inVitro, d.inVivo, d.subDephospoSites
FROM depodPhosphatase d
/*WHERE uniProtIDSub IN (%s)
LIMIT 15*/;
"""

