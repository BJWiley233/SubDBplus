USE protTest;

SELECT u.geneNamePreferred, u.uniProtID, u.proteinName, u.organism,
	u.taxid, u.geneNamesAlternative, u.alternateNames
FROM proteinsUniprot u
WHERE u.uniProtID IN ('O00187', 'O00506', 'O64517', 'P00533', 'P00736', 'P07948',
       'P08253', 'P08922', 'P09871', 'P17655', 'P18031', 'P20417',
       'P23467', 'P28482', 'P28562', 'P29350', 'P30304', 'P31749',
       'P42574', 'P42685', 'P49137', 'P51452', 'P55211', 'P60484',
       'P62136', 'P67775', 'Q12913', 'Q12923', 'Q13362', 'Q13546',
       'Q16828', 'Q80W65', 'Q96SB4', 'Q9UM73', 'Q9Y572');

SELECT * FROM proteinsUniprot WHERE uniProtID = 'P31749';