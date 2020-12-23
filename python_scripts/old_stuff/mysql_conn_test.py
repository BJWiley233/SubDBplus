import mysql.connector

cnx = mysql.connector.connect(user='root', password='**',
                              host='127.0.0.1',
                              database='merops')
cursor = cnx.cursor()   

cols = ["Protease", "Substrate_name", "Substrate_formula"]                        

query = ("SELECT %s, %s, %s  \
			 FROM merops.Substrate_search \
			 LIMIT 10;" % (cols[0], cols[1], cols[2]))


cursor.execute(query)

query = "SELECT Protease, Substrate_name, Substrate_formula FROM merops.Substrate_search LIMIT 10;"
cursor.execute(query)

query = "SELECT Protease, %s, %s FROM merops.Substrate_search LIMIT 10;"
cursor.execute(query, ("Substrate_name", "Substrate_formula"))

for (Protease, Substrate_name, Substrate_formula) in cursor:
	print("{}, {}, {}".format(Protease, Substrate_name, Substrate_formula))


cnx.close()
