SELECT  DISTINCT i.interactionID, i.interactorA, i.interactorB, i.taxidA, i.taxidB,
	i.geneNamePreferredA as nameA, i.geneNamePreferredB as nameB, i.organismA, i.organismB, 
    i.alternateNamesA as altProtNamesA, i.alternateNamesB as altProtNamesB, lower(i.direction) as direction,
    i.geneNamesAlternativeA as altGeneNamesA, i.geneNamesAlternativeB as altGeneNamesB,
    i.pubmedID, i.detectionMethod, i.interactionType, sourceDB,
    CASE
		WHEN i.interactorA REGEXP "^EBI-" THEN "Molecule"
        WHEN i.interactorA REGEXP "^[0-9]" THEN "Molecule"
        ELSE "Protein"
	END as nodeTypeA,
    CASE
		WHEN i.interactorB REGEXP "^EBI-" THEN "Molecule"
        WHEN i.interactorB REGEXP "^[0-9]" THEN "Molecule"
        ELSE "Protein"
	END as nodeTypeB
FROM protTest.intact i
WHERE interactionID IN ("EBI-762547");


SELECT interactionID, interactorA, interactorB, geneNamePreferredA, geneNamePreferredB,
	CASE
		WHEN i.interactorA REGEXP "^EBI-" THEN "Molecule"
        WHEN i.interactorA REGEXP "^[0-9]" THEN "Molecule"
        ELSE "Protein"
	END as nodeTypeA,
    CASE
		WHEN i.interactorB REGEXP "^EBI-" THEN "Molecule"
        WHEN i.interactorB REGEXP "^[0-9]" THEN "Molecule"
        ELSE "Protein"
	END as nodeTypeB
FROM intact i
WHERE length(interactorA)>6 or length(interactorB)>6;


SELECT *
FROM protTest.intact i
WHERE direction = "A->A";

SET SQL_SAFE_UPDATES = 0;
UPDATE intact
SET direction = "A->A"
WHERE interactorA = interactorB;
