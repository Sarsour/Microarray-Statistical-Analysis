# Load swirl data set
dat <- data(swirl)


# MvA plot of array 3
swirl.array3 <- swirl[,3]
print(dim(swirl.array3))

maPlot(swirl.array3, main = "Swirl Array 3 - MvA Plot", lines.func = NULL, legend.func = NULL)


# Normalize array 3 by global median location normalization
swirl.Norm.Array <- maNorm(swirl.array3, norm = c("median"))


# MvA plot of normalized array
maPlot(swirl.Norm.Array, main = "Swirl Array 3 (Normalized) - MvA Plot", lines.func = NULL, legend.func = NULL)


# MvA of normalized data with loess
swirl.Loess.Norm.Array <- maNorm(swirl.array3, norm = c("loess"))
maPlot(swirl.Loess.Norm.Array, main = "Swirl Array 3 (Normalized - Loess) - MvA Plot", lines.func = NULL, legend.func = NULL)


# Import cDNA array data
dir.path <- "C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab4"
a.cdna <- read.GenePix(path=dir.path,name.Gf = "F532 Median",name.Gb ="B532 Median", name.Rf = "F635 Median", name.Rb = "B635 Median",name.W ="Flags")


# Normalize both arrays and plot MvA for 3 different methods
patient1 <- a.cdna[ ,1]
patient1.Loess.Normalization <- maNorm(patient.one, norm = c("printTipLoess"))
patient1.MAD <- maNorm(patient.one, norm = c("scalePrintTipMAD"))
patient2 <- a.cdna[ ,2]
patient2.Loess.Normalization <- maNorm(patient.one, norm = c("printTipLoess"))
patient2.MAD <- maNorm(patient.one, norm = c("scalePrintTipMAD"))

par(mfrow = c(3,1))
maPlot(patient1, main = "Patient 1 MvA (No Normalization)", lines.func = NULL, legend.func = NULL)
maPlot(patient1.Loess.Normalization, main = "Patient 1 MvA (Print-tip Loess Normalization)", lines.func = NULL, legend.func = NULL)
maPlot(patient1.MAD, main = "Patient 1 MvA (Scale Print-tip Normalization using the MAD)", lines.func = NULL, legend.func = NULL)

par(mfrow = c(3,1))
maPlot(patient2, main = "Patient 2 MvA (No Normalization)", lines.func = NULL, legend.func = NULL)
maPlot(patient2.Loess.Normalization, main = "Patient 2 MvA (Print-tip Loess Normalization)", lines.func = NULL, legend.func = NULL)
maPlot(patient2.MAD, main = "Patient 2 MvA (Scale Print-tip Normalization using the MAD)", lines.func = NULL, legend.func = NULL)


# Data matrix output
combined.Data <- cbind(patient1, patient2)
combined.Data.Loess.Normalization <- cbind(patient1.Loess.Normalization, patient2.Loess.Normalization)
combined.Data.MAD.Normalization  <- cbind(patient1.MAD, patient2.MAD)

probe.ids <- combined.Data@maGnames@maLabels
Loess.Data <- maM(combined.Data.Loess.Normalization)
MAD.Data <- maM(combined.Data.MAD.Normalization)

dimnames(Loess.Data)[[1]] <- probe.ids
dimnames(MAD.Data)[[1]] <- probe.ids
print(dim(Loess.Data))
print(dim(MAD.Data))


# Read Affymetrix data
library(affy)
library(limma)
library(simpleaffy)
library(affyPLM)
library(fpc)

fns <- sort(list.celfiles(path = dir.path, full.names = TRUE))
data.affy <- ReadAffy(filenames = fns, phenoData = NULL)


# Normalized data matrices for affymetrix data
affy.justMAS <- justMAS(data.affy)
affy.rma <- call.exprs(data.affy, "rma")
print(dim(affy.justMAS))
print(dim(affy.rma))


# Calculate correlation between the arrays
rma.cor <- cor(exprs(affy.rma))
MAS.cor <- cor(exprs(affy.justMAS))

print(rma.cor)
print(MAS.cor)

