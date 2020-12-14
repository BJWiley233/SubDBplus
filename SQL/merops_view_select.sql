SELECT * FROM protTest.merops_select
WHERE proteaseName LIKE '%Beta-secretase 1%';

/* Original MEROPS column */
SELECT * FROM Substrate_search
WHERE Protease LIKE '%BACE1%';

/* Joined by merops_id to geneNamePreferred */
SELECT * FROM merops_select
WHERE proteaseGenePreferred LIKE '%BACE1%';