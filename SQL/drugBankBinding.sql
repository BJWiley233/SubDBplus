SELECT DISTINCT j.`chain` as 'UniProt Protein Chain', pb.*, sel.drugBankID
FROM uniprotPdbJoin j
INNER JOIN 
 pdbBindingSites pb
 ON j.pdbID = pb.pdbID
LEFT JOIN 
 (SELECT DISTINCT d.drugShort, d.drugBankID FROM pdbDrugBank d) sel
  ON pb.ligandShort = sel.drugShort
WHERE j.uniProtID = '%s'