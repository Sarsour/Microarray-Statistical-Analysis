# Downloading data set and checking out the data
dat <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\HW1\\renal_cell_carcinoma.txt", header=T, row.names=1, na.strings="NA", blank.lines.skip=F)

num.arrays <- ncol(dat)
num.probesets <- nrow(dat)

print(num.arrays)
print(num.probesets)

print(dim(dat))

dat <- as.data.frame(dat)


# Labeling the columns as desired

ann <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\HW1\\renal_carcinoma_annotation.txt", 
	na.strings="NA", blank.lines.skip=F)

gsm.id <- colnames(dat)
nt.identity <- as.character(ann[,9])
nt.identity <- sub(";", "", nt.identity)
nt.identity <- sub("[(]", "", nt.identity)
nt.identity <- sub("[)]", "", nt.identity)

print(gsm.id)
print(nt.identity)

new.header <- paste0(gsm.id, "_", nt.identity)
print(new.header)

colnames(dat) <- new.header

print(colnames(dat))


# Correlation Heat Map Plot

library(gplots)
dat.cor <- cor(dat)
 
layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10)
image(dat.cor,main="HGU133A Correlation Plot (Heat Map)",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis =0.9,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis =0.9,las=2) 

image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)


# Hierarchical Clustering Dendrogram

dat.t <- t(dat) 
dat.dist <- dist(dat.t) 
dat.clust <- hclust(dat.dist, method = "complete", members = NULL)

options(scipen=5)
plot(dat.clust, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = "Samples", ylab = "Height")


# CV vs Mean Plot

dat.mean <- apply(log2(dat),2,mean)
dat.sd <- sqrt(apply(log2(dat),2,var))
dat.cv <- dat.sd/dat.mean

plot(dat.mean,dat.cv,main="HGU133A Data Set\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="green",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)


# Average Correlation Plot
 
dat.cor <- cor(dat) 
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of HGU133A Data Set",axes=F)
points(dat.avg,bg="Blue",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")


# Remove Outliers GSM146798_Normal and GSM146799_Tumor

library(impute)

dat <- dat[, -c(10, 19)]
print(colnames(dat))


# Isolate and plot Profile Plot for AQP2 and KNG1

KNG1.1 <- as.numeric(dat["206054_at", ])
KNG1.2 <- as.numeric(dat["217512_at", ])
AQP2  <- as.numeric(dat["206672_at", ])

plot(range(1:20), range(AQP2, na.rm=TRUE), type = "n", xlab = "",
ylab = "Expression", main = "AQP2 - Gene Profile Plot", axes=F)
axis(side=1,at=c(1:20),labels=colnames(dat),cex.axis=0.4,las=2)
axis(side=2)
lines(c(1:20), AQP2, lwd=2)
grid(col = "grey")
mtext(side = 1, text = "Samples", line = 4)

plot(range(1:20), range(KNG1.1, na.rm=TRUE), type = "n", xlab = "",
ylab = "Expression", main = "KNG1 (206054_at) - Gene Profile Plot", axes=F)
axis(side=1,at=c(1:20),labels=colnames(dat),cex.axis=0.4,las=2)
axis(side=2)
lines(c(1:20), KNG1.1, lwd=2)
grid(col = "grey")
mtext(side = 1, text = "Samples", line = 4)

plot(range(1:20), range(KNG1.2, na.rm=TRUE), type = "n", xlab = "",
ylab = "Expression", main = "KNG1 (206672_at) - Gene Profile Plot", axes=F)
axis(side=1,at=c(1:20),labels=colnames(dat),cex.axis=0.4,las=2)
axis(side=2)
lines(c(1:20), KNG1.2, lwd=2)
grid(col = "grey")
mtext(side = 1, text = "Samples", line = 4)


# Assess accuracy of the missing value imputation - remove value

dat.matrix <- as.matrix(dat)
orig.value <- dat.matrix["206054_at", "GSM146784_Normal"] 
print(orig.value)

dat.matrix["206054_at", "GSM146784_Normal"] <- NA
na.value <- dat.matrix["206054_at", "GSM146784_Normal"]
print(na.value)


# Estimate missing value

dat.matrix <- impute.knn(dat.matrix, 6)
new.value <- dat.matrix$data["206054_at", "GSM146784_Normal"]
dat.knn <- dat.knn$data
print(new.value)


# Calculate relative error

absolute.error <- abs(orig.value - new.value)
rel.error <- (absolute.error/orig.value)*100
print(rel.error)


# SVD Imputation to impute missing value

library(pcaMethods)
dat.matrix2 <- as.matrix(dat) 
dat.matrix2["206054_at", "GSM146784_Normal"] <- NA
dat.matrix2 <- pca(dat.matrix2, method = "svdImpute", nPcs = 9) 
dat.pca.final <- completeObs(dat.matrix2) 
print(dat.pca.final["206054_at", "GSM146784_Normal"])
print(new.value)


# Gene profile plot of the original and new values

orig.set <- dat["206054_at",] 
knn.set <- dat.knn["206054_at",]
pca.set <- dat.pca.final["206054_at",] 
all.set <- rbind(orig.set, knn.set, pca.set) 
 
plot(c(1, ncol(dat)), range(all.set), type = "n", main = "Gene Profile Plot\nActual and Predicted Values\n(GSM146784_Normal - 206054_at)", xlab = "", ylab = "Expression Intensity", axes = F) 
axis(side = 1, at = 1:length(dat), labels = dimnames(all.set)[[2]], cex.axis = 0.4, las = 2) 
axis(side = 2)
mtext(side = 1, text = "Samples", line = 4)
for(i in 1:length(all.set)) { 
	lines(c(1:ncol(all.set)), as.numeric(all.set[i, ]), col = i+3, lwd = 2) 
} 
