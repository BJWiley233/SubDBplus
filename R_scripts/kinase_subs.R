data <- read.csv("~/JHU_Fall_2020/Biological_DBs/Project/PhosphoSitePlus_Substrates_of_Kinases/gene_attribute_matrix.txt",
                   header=T, sep='\t', stringsAsFactors = F)

data2 <- data[-1, -2]
data3 <- data2
rownames(data3) <- data3$X.
data3 <- data3[, -1]

data4 <- data3
data4[2:dim(data3)[2]] <- sapply(data4[2:dim(data3)[2]], as.numeric)

# THE ROWS ARE THE SUBSTRATES!!!!!!!!!!!!!!!!!!!!!!!
# THE COLUMNS ARE THE KINASES
# REMEMBER THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
row.names(data4)[data4$PIK3CA==1]
data4$GeneSym[data4$PIK3CA==1]

row.names(data4)[data4$CDK12==1]


"MTOR" %in% row.names(data4)[data4$AKT1==1]
"AKT1" %in% row.names(data4)[data4$PDPK1==1]
"MARK2" %in% row.names(data4)[data4$GSK3B==1]

substr.of.CSK = row.names(data4)[data4$CSK==1]
data4$GeneSym[data4$CSK==1]
data4[substr.of.CSK, "GeneSym"]

row.names(data4)[data4$PRKAB1==1]
data4$PRKAB

row.names(data4)[data4$CAMK2D==1]
row.names(data4)[data4$SRC==1]
row.names(data4)[data4$FGFR1==1]

kinase_subs <- read.table("/home/coyote/JHU_Fall_2020/Biological_DBs/Project/PhosphoSitePlus_Substrates_of_Kinases/Phosphorylation_site_dataset",
                          skip=3, sep="\t", header = T)

kinase_subs$ORGANISM
library(dplyr)
kinase_subs %>% distinct(ORGANISM)
kinase_subs[kinase_subs$ORGANISM=="fruit fly", ]
library(igraph)

intact <- read.table("/home/coyote/JHU_Fall_2020/Biological_DBs/Project/data/intact_cut_direction.txt",
                     sep="\t", )
