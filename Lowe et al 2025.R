# Rhizarian high-rank assemblage change through time in the Scotia Sea: a paleo-genomics approach (IODP Exp. 382)
# Lowe et al., 2025
# Analysis Code

#===
# Load libraries and data
#====
library(vegan)
library(adespatial)



RSD <- read.csv("RSD_RA.csv")
PR2 <- read.csv("PR2_RA.csv")

#U1534
PR2U4 <- PR2[c(1:15),]
PR2U4spe <- PR2U4[,c(8:11)]
RSDU4 <- RSD[c(1:15),]
RSDU4spe <- RSDU4[,c(8:11)]

#U1536
PR2U6 <- PR2[c(17:29),]
PR2U6spe <- PR2U6[,c(8:11)]
RSDU6 <- RSD[c(17:29),]
RSDU6spe <- RSDU6[,c(8:11)]

#U1538
PR2U8 <- PR2[c(30:61),]
PR2U8spe <- PR2U8[,c(8:11)]
RSDU8 <- RSD[c(30:61),]
RSDU8spe <- RSDU8[,c(8:11)]

# Data Transformations
PR2U4spe.hel <- decostand(PR2U4spe, "hellinger")
RSDU4spe.hel <- decostand(RSDU4spe, "hellinger")
PR2U6spe.hel <- decostand(PR2U6spe, "hellinger")
RSDU6spe.hel <- decostand(RSDU6spe, "hellinger")
PR2U8spe.hel <- decostand(PR2U8spe, "hellinger")
RSDU8spe.hel <- decostand(RSDU8spe, "hellinger")

PR2U4spe.ch <- dist.ldc(PR2U4spe, "chord")
RSDU4spe.ch <- dist.ldc(RSDU4spe, "chord")
PR2U6spe.ch <- dist.ldc(PR2U6spe, "chord")
RSDU6spe.ch <- dist.ldc(RSDU6spe, "chord")
PR2U8spe.ch <- dist.ldc(PR2U8spe, "chord")
RSDU8spe.ch <- dist.ldc(RSDU8spe, "chord")

#====
# Pearson Correlation
#====

U4correlation_Acan <- cor.test(PR2U4spe$Acantharea, RSDU4spe$Acantharea, method = 'pearson')
U4correlation_Phae <- cor.test(PR2U4spe$Phaeodarea, RSDU4spe$Phaeodarea, method = 'pearson')
U4correlation_Poly <- cor.test(PR2U4spe$Polycystines, RSDU4spe$Polycystines, method = 'pearson')
U4correlation_Taxo <- cor.test(PR2U4spe$Taxopodida, RSDU4spe$Taxopodida, method = 'pearson')

U6correlation_Acan <- cor.test(PR2U6spe$Acantharea, RSDU6spe$Acantharea, method = 'pearson')
U6correlation_Phae <- cor.test(PR2U6spe$Phaeodarea, RSDU6spe$Phaeodarea, method = 'pearson')
U6correlation_Poly <- cor.test(PR2U6spe$Polycystines, RSDU6spe$Polycystines, method = 'pearson')
U6correlation_Taxo <- cor.test(PR2U6spe$Taxopodida, RSDU6spe$Taxopodida, method = 'pearson')

U8correlation_Acan <- cor.test(PR2U8spe$Acantharea, RSDU8spe$Acantharea, method = 'pearson')
U8correlation_Phae <- cor.test(PR2U8spe$Phaeodarea, RSDU8spe$Phaeodarea, method = 'pearson')
U8correlation_Poly <- cor.test(PR2U8spe$Polycystines, RSDU8spe$Polycystines, method = 'pearson')
U8correlation_Taxo <- cor.test(PR2U8spe$Taxopodida, RSDU8spe$Taxopodida, method = 'pearson')

U4correlation_Acan
U4correlation_Phae
U4correlation_Poly
U4correlation_Taxo

U6correlation_Acan
U6correlation_Phae
U6correlation_Poly
U6correlation_Taxo

U8correlation_Acan
U8correlation_Phae
U8correlation_Poly
U8correlation_Taxo
cor
#====
# Cluster Analysis
#====
PR2U4spe.ch.UPGMA <- hclust(PR2U4spe.ch, method = "average")
RSDU4spe.ch.UPGMA <- hclust(RSDU4spe.ch, method = "average")
PR2U6spe.ch.UPGMA <- hclust(PR2U6spe.ch, method = "average")
RSDU6spe.ch.UPGMA <- hclust(RSDU6spe.ch, method = "average")
PR2U8spe.ch.UPGMA <- hclust(PR2U8spe.ch, method = "average")
RSDU8spe.ch.UPGMA <- hclust(RSDU8spe.ch, method = "average")

#====
# Changepoint Analysis
#====
