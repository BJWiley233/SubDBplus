USE protTest;

CREATE OR REPLACE VIEW merops_select AS
SELECT sel.cleavageID, sel.Ref, sel.cleavage_type, sel.Substrate_formula,
	sel.Substrate_name as substrateName, m.geneNamePreferred as substrateGenePreferred,
    m.geneNamesAlternative as subGeneAlt, m.taxid as substrateTax,
	sel.substrateUniprot, sel.substrateOrganism, 
	sel.Protease, sel.`code`, sel.proteaseUniprot, sel.geneNamePreferred as proteaseGenePreferred,
    sel.proteinName as proteaseName, sel.alternateNames as proteaseAltNames, 
    sel.taxid as proteaseTaxId, sel.organism as proteaseOrganism, sel.meropsID
FROM
(SELECT s.cleavageID, s.Ref, s.cleavage_type, s.Substrate_formula,
	s.Substrate_name, SUBSTRING(s.Uniprot, 1, 6) AS substrateUniprot, s.organism as substrateOrganism, 
	s.Protease, s.`code`, p.uniProtID as proteaseUniprot, p.geneNamePreferred, p.proteinName, 
    p.alternateNames, p.taxid, p.organism, p.meropsID
FROM protTest.Substrate_search s
INNER JOIN 
	protTest.proteinsUniprot p
    ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%', p.organism, '%')) as sel
LEFT JOIN
	meropsSubsUniprot m
    ON sel.substrateUniprot = m.uniProtID
WHERE sel.substrateUniprot IS NOT NULL;