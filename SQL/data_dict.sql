SELECT
tb.COLUMN_NAME AS 'Field Name',
tb.COLUMN_TYPE AS 'Data Type',
tb.IS_NULLABLE AS 'Allow Empty',
tb.EXTRA AS 'PK',
tb.COLUMN_COMMENT AS 'Field Description',
tb.table_schema,
tb.table_name
FROM
`INFORMATION_SCHEMA`.`COLUMNS` as tb
WHERE table_schema = 'protTest';


SELECT
tb.COLUMN_NAME AS 'Element name',
tb.table_name AS 'MySQL table name',
'' AS Description, '' AS Source,
tb.COLUMN_TYPE AS 'Data type',
'' AS 'Neo4j node or edge',
'' AS 'Neo4j element name',
'' AS 'Neo4j data type'
FROM
`INFORMATION_SCHEMA`.`COLUMNS` as tb
WHERE table_schema = 'protTest'
ORDER BY LOWER(tb.table_name) ASC;