# SubDBplus
## Updated

### Abstract
An important question to ask in drug discovery when targeting a gene or expression of that gene is
what effect can my drug have on that protein’s downstream and upstream pathways. But how do we
start to analyze various downstream pathways if we don’t have the entire picture of all the substrates,
either proteins, DNA, RNA or small molecules that our protein has or interacts with. More specifically
it would be important to know where in our protein, i.e. which domains or residues are specific to each
interaction between protein and substrate. I have developed a database that as its primary function will
give the researcher, whether they be a molecular biologist, computational chemist, or in lead-to-target
drug discovery, a list of substrates for a particular protein of interest. That database, called SubDBplus
(SubDB+) is implemented in both a relational MySQL data warehouse as well as the graphing NoSQL
database Neo4j so that in addition to substrates, up and downstream pathways can also be queried. The
database forms a collaboration between 4 different substrate or binary interactions databases including
phosphorylation and dephosphorylation databases, proteolysis, and a non-specific or various interaction
types database that includes but not limited to ubiquitinases, transferases, methyltransferases,
acetylases, and G-protein coupled receptors. Up until writing this paper, I was not aware of any
interfaces that included a simple way to search a protein’s substrates and include a graphing model.
However during the development of the database and writing the manuscript, I have also discovered a
similar database specific to phosphorylation, RegPhos, which is also implemented in MySQL but to my
knowledge does not have a graphing layer to the caliber of Neo4j. During the development of SubDB+
and querying of the interface I have learned that ideas for improvements and extended releases of a
database come from using the database, seeing where things can be better as well as seeing limitations
and where certain functionalities are best supported in different application layers.

