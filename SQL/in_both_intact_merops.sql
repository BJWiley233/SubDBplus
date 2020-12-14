USE protTest;

/* Also gets substrate Gene name preferred if in Uniprot table */
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