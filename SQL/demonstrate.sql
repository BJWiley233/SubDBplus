/* Which kinase has most unique substrates */
SELECT * 
FROM kinasePhosphoSitePlus
WHERE uniProtIDKin = 
	(SELECT sel1.uniProtIDKin
     FROM (SELECT DISTINCT uniProtIDKin, uniProtIDSub
							FROM kinasePhosphoSitePlus) as sel1
					  GROUP BY uniProtIDKin
					  ORDER BY COUNT(uniProtIDKin) DESC
					  LIMIT 1);

SELECT * 
FROM uniprotPdbJoin
WHERE uniProtID = 'P17612';

SELECT * 
FROM intact
WHERE interactorA = 'P17612';

/* Which peptidase has most unique substrates */
SELECT proteaseUniprot, proteaseGenePreferred, count(proteaseUniprot)
/* merops_select is a complex view */
FROM protTest.merops_select
GROUP BY proteaseUniprot
ORDER BY count(proteaseUniprot) DESC;

SELECT * FROM merops_select m
WHERE m.proteaseName LIKE '%Beta-secretase 1%';

SELECT * FROM protTest.merops_select m
WHERE m.proteaseGenePreferred LIKE '%BACE1%';




