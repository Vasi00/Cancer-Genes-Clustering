# libraries
library(dplyr)
library(caret)
library(psych)

# Reading the data and viewing it
Cancer_data<-read.table("Cancer data.txt",sep="",row.names = 1, na.strings = c("Null"))
View(Cancer_data)

# Structure and summary of dataset
str(Cancer_data)
summary(Cancer_data)
describe(Cancer_data)

# Setting row names according to genes (copy names from row-file)
# row names represent different gene names
data1<-as.matrix(Cancer_data)
row.names(data1) <- c("18s","18S","18S","18s","ABCA1","ABCA10","ABCA12","ABCA13","ABCA2","ABCA3","ABCA4","ABCA5","ABCA6","ABCA7","ABCA8",	"ABCA9",	"ABCB1",	"ABCB10",	"ABCB11",	"ABCB4",	"ABCB5",	"ABCB6",	"ABCB7",	"ABCB8",	"ABCB9","ABCC1",	"ABCC10",	"ABCC11",	"ABCC12",	"ABCC2",	"ABCC3",	"ABCC4",	"ABCC5",	"ABCC6",	"ABCC8",	"ABCC9","ABCD1",	"ABCD2",	"ABCD3",	"ABCD4",	"ABCE1",	"ABCF1",	"ABCF2",	"ABCF3",	"ABCG1",
                      "ABCG2",	"ABCG4", "18s",	"18S",	"18S",	"18s",	"ABCA1",	"ABCA10",	"ABCA12",	"ABCA13",	"ABCA2",	"ABCA3",	"ABCA4",	"ABCA5","ABCA6",	"ABCA7",	"ABCA8",	"ABCA9",	"ABCB1",	"ABCB10",	"ABCB11",	"ABCB4",	"ABCB5",	"ABCB6",	"ABCB7",	"ABCB8",	"ABCB9",	"ABCC1",	"ABCC10",	"ABCC11",	"ABCC12",	"ABCC2",	"ABCC3",	"ABCC4",	"ABCC5",	"ABCC6",	"ABCC8",	"ABCC9","ABCD1",	"ABCD2",	"ABCD3",	"ABCD4",
                      "ABCE1",	"ABCF1",	"ABCF2",	"ABCF3",	"ABCG1",	"ABCG2",	"ABCG4",	"ABCG5","ABCG8",	"ABL1",	"ACTB",	"AHR",	"AKR1C1 AKR1C2",	"AKT1",	"ANXA1","ANXA4",	"APAF1",	"APC",	"APEX1",	"APEX2", "APOE",	"AQP7",	"AQP9",	"ASAH1",	"ASAH2",	"ASAH3",	"ATM",	"ATOX1",	"ATP1A1",	"ATP1B1",	"ATP6V0C",	"ATP7A", "ATP7B",	"ATP8B1",	"ATR",	"AURKA",	"BAD",	"BAG1",	"BAG3",	"BAG4",	"BAK1",	"BAX",	"BCL2",	"BCL2A1",	"BCL2L1",
                      "BCR",	"BID",	"BIRC2",	"BIRC3",	"BIRC4",	"BIRC5",	"BIRC6",	"BNIP3",	"BNIP3L",	"BRCA1",	"BRCA2",	"C8orf33","CASP3",	"CASP8",	"CASP9",	"CCL2",	"CCND1",	"CCNE1",	"CCNH",	"CCT8",	"CDC40",	"CDC42",	"CDH1",	"CDK2",	"CDK4", "CDK7",	"CDKN1A",	"CDKN1B",	"CDKN2A",	"CFL1",	"CFLAR",	"CFTR",	"CHEK1",	"CHEK2",	"CHUK",	"CIAPIN1",	"CLDN1",	
                      "CLDN16",	"CLDN2",	"CLDN3",	"CLDN4",	"CLDN5",	"CLDN7",	"CLPTM1L",	"CLU",	"COX7A2",	"CROP",	"CTNNA1",	"CTNNB1",
                      "CYP1A2",	"CYP2A13 CYP2A6 CYP2A7",	"CYP2B6",	"CYP2C19 CYP2C8",	"CYP2C8",	"CYP2C9",	"CYP2D6",	"CYP2E1",	"CYP3A4",	
                      "CYP3A5",	"DAPK1",	"DHFR",	"DIABLO",	"E2F1",	"EGFR",	"HAGH",	"HIF1A",	"HSF1",	"HSP90AA1",	"HSPA5",	"HSPB1",	"HSPD1",
                      "HSPE1",	"HSPH1",	"HUS1",	"IGF1R",	"IL6",	"INSR",	"ITGAE",	"ITGB1",	"JUN",	"KCNMA1",	"KIT",	"KLF1",	"LAMP1",
                      "LAMP2",	"LBR",	"LIG4",	"MAP2K1",	"MAPK1",	"MAPK3",	"MAPK8",	"MBD4",	"MCL1",	"MDM2",	"MGMT",	"MKI67",	"MLH1",	
                      "MLH3",	"MMP2",	"MMP9",	"MNAT1",	"MRE11A",	"MSH2",	"MSH3",	"MSH6",	"MT1A",	"MT1B",	"MT1F",	"MT1H",	"MT1X",	"MT2A",	
                      "MT3",	"MT4",	"MTMR11",	"MUTYH",	"MVP",	"MYC",	"NFKB1",	"NFKBIA",	"NOLA2",	"NR1H2",	"NR1H3",	"NR1H4",	"NR1I2",
                      "NR1I3",	"NRAS",	"NTRK2",	"OCLN",	"OGG1",	"OVCA2",	"PARP1",	"PARP2",	"PDCD8",	"PDGFRB",	"PDK1",	"PIK3CA",	"POLB",	
                      "POLH",	"POLI",	"POLK",	"PSMA3",	"PTEN",	"RAD1",	"RAD17",	"RAD18",	"RAD23A",	"RAD23B",	"RAD50",	"RAD51",	"RAF1",
                      "RB1",	"RELA",	"RHOD",	"RPL13A",	"RPL36",	"RPL41",	"RXRB",	"S100A10",	"SEPX1",	"SFN",	"SGPP1",	"SIRT1",	"SIRT2",
                      "SIRT3",	"SIRT4",	"SIRT5",	"SIRT6",	"SIRT7",	"SLC10A1",	"SLC15A1",	"SLC15A2",	"SLC16A2",	"SLC16A3",	"SLC19A1",
                      "SLC19A2",	"SLC19A3",	"SLC1A4",	"SLC1A5",	"SLC22A1",	"SLC22A2",	"SLC25A15",	"SLC25A30",	"SLC25A5",	"SLC28A1",
                      "SLC28A3",	"SLC29A1",	"SLC29A2",	"SLC2A5",	"SLC31A1",	"SLC34A2",	"SLC3A1",	"SLC3A2",	"SLC5A4",	"SLC5A6",	"SLC7A1",
                      "SLC7A10",	"SLC7A11",	"SLC7A2",	"SLC7A3",	"SLC7A5",	"SLC7A8",	"SLC7A9",	"SLC9A3R2",	"SLCO1B3",	"SLCO4A1",	"SP1","SRC","STARD4",	"STAT1",	"STAT3",	"STAT5A",	"STAT5B",	"TAP1",	"TAP2",	
                      "TCEAL4",	"TGFA",	"TGFB1",	"TIMP1",	"TMEM109","TNF",
                      "TNFRSF10A",	"TNFSF10",	"TOP1",	"TOP2A", "TOP2B","TP53")
View(data1)

# Setting coloumn names (from col-names file)
# Column names represent Cancer cell lines/tissue
df<-as.data.frame(data1)
colnames(df)<-c("GSM839031","GSM839032","GSM839033","GSM839034","GSM839035","GSM839036","GSM839037","GSM839038","GSM839039","GSM839040","GSM839041","GSM839042","GSM839043","GSM839044","GSM839045",
                "GSM839046","GSM839047","GSM839048","GSM839049","GSM839050","GSM839051","GSM839052","GSM839053","GSM839054","GSM839055","GSM839056","GSM839057","GSM839058","GSM839059","GSM839060",
                "GSM839061","GSM839062","GSM839063","GSM839064","GSM839065","GSM839066","GSM839067","GSM839068","GSM839069")
View(df)

# missing values
table(complete.cases(Cancer_data))  # counting all values that are filled(complete) as TRUE and NA as FALSE
table(is.na(Cancer_data))           # counting total NA values

# We can't impute Null/NA values by regular mean, median, knn imputation since genes expression values are different for different genes and have no relation with other gene expression (RT-PCR) values 
# SO I will remove all the NA values
df<-na.omit(df)
View(df)

# installing preprocessCore library for preprocessing and normalization
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("preprocessCore")

library(preprocessCore)

# Quantile Normalization to make two distributions identical in statistical properties
data_mat<-data.matrix(df)
data_norm <- normalize.quantiles(data_mat,copy=TRUE)

colnames(data_norm) <- colnames(data_mat)

View(data_norm)

# Since after removing NA values, row names may have been mismatched, so we are setting them again
## Setting row names (copied from row-file)
row.names(data_norm)<- c("18S",	"18S",	"18S",	"18S",	"ABCA1",	"ABCA10",	"ABCA12",	"ABCA13",	"ABCA2",	"ABCA3",	"ABCA4",	"ABCA5",	"ABCA6",	"ABCA7",	"ABCA8",
                         "ABCA9",	"ABCB1",	"ABCB10",	"ABCB11",	"ABCB4",	"ABCB5",	"ABCB6",	"ABCB7",	"ABCB8",	"ABCB9",	"ABCC1",	"ABCC10",	"ABCC11",	"ABCC12",	
                         "ABCC2",	"ABCC3",	"ABCC4",	"ABCC5",	"ABCC6",	"ABCC8",	"ABCC9",	"ABCD1",	"ABCD2",	"ABCD3",	"ABCD4",	"ABCE1",	"ABCF1",	"ABCF2",
                         "ABCF3",	"ABCG1",	"ABCG2",	"ABCG4",	"ABCG5",	"ABCG8",	"ABL1",	"ACTB",	"AHR",	"AKR1C1 AKR1C2",	"AKT1",	"ANXA1",	"ANXA4",	"APAF1",
                         "APC",	"APEX1",	"APEX2",	"APOE",	"AQP7",	"AQP9",	"ASAH1",	"ASAH2",	"ASAH3",	"ATM",	"ATOX1",	"ATP1A1",	"ATP1B1",	"ATP6V0C",	"ATP7A",
                         "ATP7B",	"ATP8B1",	"ATR",	"AURKA",	"BAD",	"BAG1",	"BAG3",	"BAG4",	"BAK1",	"BAX",	"BCL2","BCL2A1",	"BCL2L1",	"BCR",	
                         "BID",	"BIRC2",	"BIRC3",	"BIRC4",	"BIRC5",	"BIRC6",	"BNIP3",	"BNIP3L",	"BRCA1",	"BRCA2",	"C8orf33",	
                         "CASP3",	"CASP8",	"CASP9",	"CCL2",	"CCND1",	"CCNE1",	"CCNH",	"CCT8",	"CDC40",	"CDC42",	"CDH1",	"CDK2",	"CDK4",	"CDK7",	"CDKN1A",	
                         "CDKN1B",	"CDKN2A",	"CFL1",	"CFLAR",	"CFTR",	"CHEK1",	"CHEK2",	"CHUK",	"CIAPIN1",	"CLDN1",	"CLDN16",	"CLDN2",	"CLDN3",	"CLDN4",	
                         "CLDN5",	"CLDN7",	"CLPTM1L",	"CLU",	"COX7A2",	"CROP",	"CTNNA1",	"CTNNB1",	"CYP1A2",	"CYP2A13 CYP2A6 CYP2A7",	"CYP2B6",	"CYP2C19 CYP2C8",	
                         "CYP2C8",	"CYP2C9",	"CYP2D6",	"CYP2E1",	"CYP3A4",	"CYP3A5",	"DAPK1",	"DHFR",	"DIABLO",	"E2F1",	"EGFR",	"EHBP1",	"ERBB2",	"ERCC1",	"ERCC2",	"ERCC3",	"ERCC4",	
                         "ERCC5",	"ERCC6",	"ERCC8",	"ETS1",	"F3",	"FADD",	"FAS",	"FASLG",	"FKBP1A",	"FN1",	"FOS",	"FXYD2",	"FZD1",	"GART",	"GBP1",	"GGT1",
                         "GJA1",	"GLO1",	"GPR177",	"GPX1",	"GPX2",	"GPX3",	"GPX4",	"GSK3B",	"GSR",	"GSS",	"GSTA1",	"GSTA2",	"GSTA3",	"GSTA4",	"GSTA5",
                         "GSTK1",	"GSTM2",	"GSTM3",	"GSTM4",	"GSTM5",	"GSTO1",	"GSTP1",	"GSTT1",	"GSTT2",	"GSTZ1",	"HAGH",	"HIF1A",	"HSF1",	"HSP90AA1",
                         "HSPA5",	"HSPB1",	"HSPD1",	"HSPE1",	"HSPH1",	"HUS1",	"IGF1R",	"IL6",	"INSR",	"ITGAE",	"ITGB1",	"JUN",	"KCNMA1",	"KIT",	"KLF1",
                         "LAMP1",	"LAMP2",	"LBR",	"LIG4",	"MAP2K1",	"MAPK1",	"MAPK3",	"MAPK8",	"MBD4",	"MCL1",	"MDM2",	"MGMT",	"MKI67",	"MLH1",	"MLH3",	"MMP2",
                         "MMP9",	"MNAT1",	"MRE11A",	"MSH2",	"MSH3",	"MSH6",	"MT1A",	"MT1B",	"MT1F",	"MT1H",	"MT1X",	"MT2A",	"MT3",	"MT4",	"MTMR11",	"MUTYH",
                         "MVP",	"MYC",	"NFKB1",	"NFKBIA",	"NOLA2",	"NR1H2",	"NR1H3",	"NR1H4",	"NR1I2")



# HEIRARCHIAL CLUSTERING

# Calculating distance matrix
dis <- dist(data_norm)
View(as.matrix(dis))

hc <- hclust(dis)
hc
par(cex = 0.3)
plot(hc, main = "Heirarchial Clustering", xlab = "distance", ylab = "height", font = 5, hang = -1, )

# K-MEANS CLUSTERING

df_1<- as.data.frame(data_norm)
View(df_1)

library(factoextra)  
## For optimal number of clusters using Elbow method i.e. within sum of squares (wss) and silhouette method
fviz_nbclust(df_1, kmeans, method = "wss")
fviz_nbclust(df_1, kmeans, method = "silhouette")

## Final cluster using k-means algorithm
km <- kmeans(df_1,3,nstart = 25)
km

## Visualization
kmvis <- fviz_cluster(km,df_1)
kmvis

## Purity of Clusters

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters) * 100
}
n = 256    # number of values we have in dataset after removing na values
classes = sample(39, n, replace=T)
clusters = sample(3, n, replace=T)
ClusterPurity(clusters, classes)

# Cluster purity is 6.64

