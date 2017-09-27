# Import data to table
table <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab2\\spellman.txt", header=T, row.names=1)	
print(ncol(table))
print(nrow(table))


# Isolate samples 23-46
isotable <- table[, 23:46]
print(ncol(isotable))
print(nrow(isotable))


# Correlation matrix between time points
library(gplots)
isotable.cor <- cor(isotable, use = "pairwise.complete.obs", method = "pearson")

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(isotable.cor,na.rm=T),max(isotable.cor,na.rm=T),length=10)
image(isotable.cor,main="CDC15 Experiment - Pearson Correlation Plot",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(isotable.cor)),label=dimnames(isotable.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(isotable.cor)),label=dimnames(isotable.cor)[[2]],cex.axis=0.9,las=2)

image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)


# Select a gene and impute missing values
YAL002W.vector <- as.numeric(as.vector(isotable[2,]))
YAL002W.mean <- mean(YAL002W.vector, na.rm = TRUE)
print(YAL002W.mean)

YAL002W.vector.impute <- ifelse(is.na(YAL002W.vector), mean(YAL002W.vector, na.rm = TRUE), YAL002W.vector)
print(YAL002W.vector.impute)


# Gene profile plot
dat <- read.table("C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab2\\spellman.txt",header=T,row.names=1)
dat <- as.data.frame(dat);

iso <- dat[, 23:46]

YAL002W <- iso["YAL002W", ]
plot(range(1:24), range(YAL002W, na.rm=TRUE), type = "n",  xlab = "Timepoints", ylab = "Expression",   main = "YAL002W - Gene Profile Plot", axes=F)
axis(side=1,at=c(1:24),labels=dimnames(YAL002W)[[2]],cex.axis=0.4,las=2)
axis(side=2)
lines(c(1:24), as.numeric(YAL002W), lwd=2)
grid(col = "grey")

