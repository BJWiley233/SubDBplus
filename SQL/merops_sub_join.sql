/* Will only get peptidases from UniProt tables but substrates from Merops */
SELECT s.Substrate_name, s.Uniprot AS substrate_Uniprot, s.organism as substrate_Organism, 
	s.protease, s.`code`, p.uniProtID as protease_Uniprot, p.geneNamePreferred, p.proteinName, 
    p.alternateNames, p.taxid, p.organism
FROM protTest.Substrate_search s
INNER JOIN 
	protTest.proteinsUniprot p
    ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%', p.organism, '%');
    
/* Also gets substrate Gene name preferred if in Uniprot table */
SELECT sel.Substrate_name as substrateName, prot.geneNamePreferred as substrateGenePreferred,
	sel.substrateUniprot, sel.substrateOrganism, 
	sel.protease, sel.`code`, sel.proteaseUniprot, sel.geneNamePreferred as proteaseGenePreferred,
    sel.proteinName as proteaseName, sel.alternateNames as proteaseAltNames, 
    sel.taxid as proteaseTaxId, sel.organism as proteaseOrganism
FROM
(SELECT s.Substrate_name, s.Uniprot AS substrateUniprot, s.organism as substrateOrganism, 
	s.protease, s.`code`, p.uniProtID as proteaseUniprot, p.geneNamePreferred, p.proteinName, 
    p.alternateNames, p.taxid, p.organism
FROM protTest.Substrate_search s
INNER JOIN 
	protTest.proteinsUniprot p
    ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%', p.organism, '%')) as sel
LEFT JOIN
	protTest.proteinsUniprot prot
    ON sel.substrateUniprot = prot.uniProtID;



