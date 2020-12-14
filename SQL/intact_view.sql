USE protTest;

CREATE OR REPLACE VIEW intact_select AS
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
FROM protTest.intact i;


