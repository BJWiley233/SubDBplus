USE protTest;

/* Also gets substrate Gene name preferred if in Uniprot table */
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
    ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%', p.organism, '%')) as sel
LEFT JOIN
	meropsSubsUniprot m
    ON sel.substrateUniprot = m.uniProtID
WHERE sel.cleavageID IN ('CLE0000094', 'CLE0000095', 'CLE0000099', 'CLE0000100',
       'CLE0000566', 'CLE0001516', 'CLE0001517', 'CLE0002200',
       'CLE0003178', 'CLE0084358');









/* Also gets substrate Gene name preferred if in Uniprot table */
SELECT count(sel.substrateUniprot), sel.substrateUniprot
FROM
(SELECT s.cleavageID, s.Ref, s.cleavage_type, s.Substrate_formula,
	s.Substrate_name, SUBSTRING(s.Uniprot, 1, 6) AS substrateUniprot, s.organism as substrateOrganism, 
	s.protease, s.`code`, p.uniProtID as proteaseUniprot, p.geneNamePreferred, p.proteinName, 
    p.alternateNames, p.taxid, p.organism, p.meropsID
FROM protTest.Substrate_search s
INNER JOIN 
	protTest.proteinsUniprot p
    ON p.meropsID = s.`code` AND s.organism LIKE CONCAT('%', p.organism, '%')) as sel
LEFT JOIN
	meropsSubsUniprot m
    ON sel.substrateUniprot = m.uniProtID
WHERE m.geneNamePreferred is null
GROUP BY sel.substrateUniprot;



