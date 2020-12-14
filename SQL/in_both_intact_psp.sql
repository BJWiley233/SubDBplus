USE protTest;

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
LIMIT 5)

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
LIMIT 5);