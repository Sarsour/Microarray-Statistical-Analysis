# Import data into dataframe
dat <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab3\\eisen.txt", 
	header=T, row.names=1, na.strings="NA", blank.lines.skip=F)
print(ncol(dat))
print(nrow(dat))
dat <- as.data.frame(dat)


# Read in classes
ann <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab3\\eisenClasses.txt", 
	header=T)
print(ncol(ann))
print(nrow(ann))


# Subset the data
print(dimnames(dat)[[2]])
cl <- as.character(ann[,2])
dat <- dat[,cl]
print(dimnames(dat)[[2]])
class1 <- cl[1:19]
class2 <- cl[20:39]
print(colnames(class1))
print(colnames(class2))


# Boxplot and histogram
x <- as.numeric(dat[8000,class1])
y <- as.numeric(dat[8000,class2])

x <- x[!is.na(x)]
y <- y[!is.na(y)]

xy.list <- list(x, y)

boxplot(xy.list,col=c('red','blue'),main='Gene 8000 from DLBCL cDNA 2-channel dataset',axes=F,ylab='log2(ratio intensity)')
axis(2)
axis(1,at=c(1,2),c('GC','ACT'))

par(mfrow=c(2, 1))
hist(x, col = 'red', labels = T, main = 'Gene 8000 from the DLBCL cDNA 2-channel dataset - GerminalCentre Class', axes = T, 
	ylab = 'Frequency', xlab = 'log2(ratio intensity) ranges')
hist(y, col = 'blue', labels = T, main = 'Gene 8000 from the DLBCL cDNA 2-channel dataset - Activated Class', axes = T, 
	ylab = 'Frequency', xlab = 'log2(ratio intensity) ranges')


# Calculate pooled variance
nx <- length(x)
ny <- length(y)
pool.var <- (((nx-1)*var(x)) + ((ny-1)*var(y)))/(nx+ny-2)
print(pool.var)

dif.fold <- log2(1.5)/sqrt(pool.var)
print(dif.fold)

pl.ss3 <- power.t.test(d = dif.fold, sig.level = .01, power = 0.8, type = "two.sample")
print(pl.ss3)


# Find necessary sample size (99% confidence and 80% power)
dif <- abs(mean(x)-mean(y))/sqrt(pool.var)
print(dif)

pl.ss3 <- pwr.t.test(d = dif.3fold, sig.level = .01, power = 0.8, type="two.sample")
print(pl.ss3)


# Standard Deviation
library(ssize)
library(gdata) 
data(exp.sd) 

exp.sd = apply(dat, 1, sd, na.rm = T)

hist(exp.sd, col = "Red",labels = T, main = "Standard Deviation Histogram of Genes in the Matrix", 
	axes = T, ylab = "Frequency", xlab = "Standard Deviation Ranges (for data on the log2 scale)")


# Plot gene vs sample size proportion (log2 transform)
fold.change = 3.0
sig.level = 0.05 
power = 0.8

all.size <- ssize(sd=exp.sd, delta=log2(fold.change), sig.level=sig.level, power=power) 
ssize.plot(all.size, lwd=2, col="blue", xlim=c(1,20)) 
xmax <- par("usr")[2]-1; 
ymin <- par("usr")[3] + 0.05 
title("Sample Size to Detect 3-Fold Change")
legend(x=xmax, y=ymin, legend= strsplit( paste("fold change=",fold.change,",", "alpha=", sig.level, ",", "power=",power,",", "# genes=", length(exp.sd), sep=''), "," )[[1]], xjust=1, yjust=0, cex=1.0) 


