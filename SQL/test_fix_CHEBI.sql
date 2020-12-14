/* Test fixing this in neo4j, may show fix for longer interactorA
   Short story is there is delimeter ':' that we fixed and CHEBI
   Molecules have 5 or 6 numbers after 'CHEBI:' */
SELECT DISTINCT interactorB, concat('CHEBI:', interactorB) as fixedInteractorB, nameB
FROM protTest.intact_select
WHERE interactorB REGEXP "^[0-9]{5,6}$";

/* Update it for B only leaving A still incorrect */
UPDATE intact
SET interactorB = concat('CHEBI:', interactorB)
WHERE interactorB REGEXP "^[0-9]{5,6}$";
/* Confirm */
SELECT * FROM intact WHERE interactorB REGEXP "CHEBI.*";

/* All CHEBI molecules */
SELECT DISTINCT interactorA, concat('CHEBI:', interactorA) as fixedInteractor, nameA as name_
FROM protTest.intact_select
WHERE interactorA REGEXP "^[0-9]{5,6}$";
UNION
SELECT DISTINCT interactorB, concat('CHEBI:', interactorB) as fixedInteractor, nameB as name_
FROM protTest.intact_select
WHERE interactorB REGEXP "^[0-9]{5,6}$"

