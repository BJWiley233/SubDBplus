USE protTest;

(SELECT DISTINCT u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, d.uniProtIDPPase
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorA
JOIN 
	( SELECT DISTINCT d.uniProtIDPPase
      FROM depodPhosphatase d 
	) as d
    ON u.uniProtID = d.uniProtIDPPase
ORDER BY RAND()
LIMIT 5) 

UNION ALL

(SELECT DISTINCT u.geneNamePreferred, u.uniProtID, u.headProteinFamily, i.interactorA,
	i.geneNamePreferredA, i.interactorB, i.geneNamePreferredB, i.interactionType, 
    i.direction, i.interactionID, d.uniProtIDPPase
FROM proteinsUniprot u
INNER JOIN intact i 
	ON u.uniProtID = i.interactorB
JOIN 
	( SELECT DISTINCT d.uniProtIDPPase
      FROM depodPhosphatase d 
	) as d
    ON u.uniProtID = d.uniProtIDPPase
ORDER BY RAND()
LIMIT 5);