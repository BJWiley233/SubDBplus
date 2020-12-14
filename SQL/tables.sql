DROP DATABASE IF EXISTS protTest;
CREATE DATABASE IF NOT EXISTS protTest;
USE protTest;

/* maybe call it uniprotProteins */
CREATE TABLE IF NOT EXISTS proteins  (
	uniProtID CHAR(15) NOT NULL PRIMARY KEY, 
    entryName VARCHAR(45) NOT NULL, 
    reviewed BOOL NOT NULL,
    proteinName VARCHAR(255) NOT NULL,
    alternateNames JSON, /*test*/
    geneNamePreferred VARCHAR(45) NOT NULL,
    geneNamesAlternative JSON, 
    geneNamesAll JSON,
    interactsWith JSON,
    organismID INT,
    organism VARCHAR(45),
	proteinFamilies VARCHAR(255),
    length INT,
    headProteinFamily VARCHAR(45) NOT NULL
);

/*
INSERT INTO proteins 
	VALUES(
		"Q63088885", 
        "ANM1_MOUSE",
        TRUE,
        "Protein arginine N-methyltransferase 1",
        NULL,
        "Prmt1",
		NULL, 
		NULL,
		NULL,
		9606,
		"Homo sapiens (Human)",
		"Methyltransferase superfamily, RRP8 family",
		456,
		"methyltransferase"
);
*/

/* Only for Intact */
CREATE TABLE IF NOT EXISTS intactSubstrates  (
	uniProtID CHAR(15) NOT NULL PRIMARY KEY, 
    entryName VARCHAR(45) NOT NULL, 
    reviewed BOOL NOT NULL,
    proteinName VARCHAR(255) NOT NULL,
    alternateNames JSON, /*test*/
    geneNamePreferred VARCHAR(45) NOT NULL,
    geneNamesAlternative JSON, 
    geneNamesAll JSON,
    interactsWith JSON,
    organismID INT,
    organism VARCHAR(45),
	proteinFamilies VARCHAR(255),
    length INT,
    headProteinFamily VARCHAR(45) NOT NULL
);
