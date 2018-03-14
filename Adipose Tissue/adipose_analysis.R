# Load GenePix files
library(marray)
dir.path <- "C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\HW2"
dat <- read.GenePix(path=dir.path,name.Gf = "F532 Median",name.Gb ="B532 Median", name.Rf = "F635 Median", name.Rb = "B635 Median",name.W ="Flags")
print(dim(dat))


# Normalize each array using these methods: median global, loess, print-tip-group-loess and plot MvA
dat.array1 <- dat[,1]
array1.Median.Global <- maNorm(dat.array1, norm = c("median"))
array1.Loess <- maNorm(dat.array1, norm = c("loess"))
array1.PrintTip <- maNorm(dat.array1, norm = c("printTipLoess"))
par(mfrow = c(4,1))
maPlot(dat.array1, main = "Array 1 - MvA - No Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array1.Median.Global, main = "Array 1 - MvA - Median Global Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array1.Loess, main = "Array 1 - MvA - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array1.PrintTip, main = "Array 1 - MvA - Print Tip Group Normalization", lines.func = NULL, legend.func = NULL)

dat.array2 <- dat[,2]
array2.Median.Global <- maNorm(dat.array2, norm = c("median"))
array2.Loess <- maNorm(dat.array2, norm = c("loess"))
array2.PrintTip <- maNorm(dat.array2, norm = c("printTipLoess"))
par(mfrow = c(4,1))
maPlot(dat.array2, main = "Array 2 - MvA - No Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array2.Median.Global, main = "Array 2 - MvA - Median Global Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array2.Loess, main = "Array 2 - MvA - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array2.PrintTip, main = "Array 2 - MvA - Print Tip Group Normalization", lines.func = NULL, legend.func = NULL)

dat.array3 <- dat[,3]
array3.Median.Global <- maNorm(dat.array3, norm = c("median"))
array3.Loess <- maNorm(dat.array3, norm = c("loess"))
array3.PrintTip <- maNorm(dat.array3, norm = c("printTipLoess"))
par(mfrow = c(4,1))
maPlot(dat.array3, main = "Array 3 - MvA - No Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array3.Median.Global, main = "Array 3 - MvA - Median Global Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array3.Loess, main = "Array 3 - MvA - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array3.PrintTip, main = "Array 3 - MvA - Print Tip Group Normalization", lines.func = NULL, legend.func = NULL)

dat.array4 <- dat[,4]
array4.Median.Global <- maNorm(dat.array4, norm = c("median"))
array4.Loess <- maNorm(dat.array4, norm = c("loess"))
array4.PrintTip <- maNorm(dat.array4, norm = c("printTipLoess"))
par(mfrow = c(4,1))
maPlot(dat.array4, main = "Array 4 - MvA - No Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array4.Median.Global, main = "Array 4 - MvA - Median Global Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array4.Loess, main = "Array 4 - MvA - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(array4.PrintTip, main = "Array 4 - MvA - Print Tip Group Normalization", lines.func = NULL, legend.func = NULL)


# Density plots of log ratios for each normalization method (array 4)
new.dat.array4 = na.omit(maM(dat.array4))
new.array4.Median.Global = na.omit(maM(array4.Median.Global))
new.array4.Loess = na.omit(maM(array4.Loess))
new.array4.PrintTip = na.omit(maM(array4.PrintTip))

plot(density(new.dat.array4), main = "Array 4 - Density Plot", ylim=c(0, 1),xlim=c(-7, 10),col="red")
lines(density(new.array4.Median.Global), col = "blue")
lines(density(new.array4.Loess), col = "purple") 
lines(density(new.array4.PrintTip), col = "green") 

leg.txt <- c("No normalization", "Median Global Normalization", "Loess Normalization", "Print-Tip-Group Normalization")
legend(1, 0.5, legend = leg.txt, lty = c(1, 1), lwd = c(2.5, 2.5),col = c("red", "blue", "purple", "green"))


# Calculate correlation after CY5 extraction, and calculate global median normalization

array1.diff <- log2((maRf(dat[ , 1])) - (maRb(dat[ , 1])))
array2.diff <- log2((maRf(dat[ , 2])) - (maRb(dat[ , 2])))
array3.diff <- log2((maRf(dat[ , 3])) - (maRb(dat[ , 3])))
array4.diff <- log2((maRf(dat[ , 4])) - (maRb(dat[ , 4])))

arrays.con <- cbind(array1.diff, array2.diff, array3.diff, array4.diff)
arrays.med <- apply(arrays.con, 2, median, na.rm = T)
arrays.norm <- sweep(arrays.con, 2, arrays.med)
colnames(arrays.norm) <- c("Array 1", "Array 2", "Array 3", "Array 4")

median(arrays.norm[ , 1], na.rm = T) 
median(arrays.norm[ , 2], na.rm = T)
median(arrays.norm[ , 3], na.rm = T)
median(arrays.norm[ , 4], na.rm = T)


# Spearman's rank correlation and plot scatter plot matrix for the two normalizations
samp.1 <- maM(array1.Loess)
samp.2 <- maM(array2.Loess)
samp.3 <- maM(array3.Loess)
samp.4 <- maM(array4.Loess)
names <- c("Array 1", "Array 2", "Array 3", "Array 4")
data.loess <- cbind(samp.1, samp.2, samp.3, samp.4)
colnames(data.loess) <- make.names(names, unique = TRUE)
data.loess.sp <- cor(data.loess, use = "complete.obs", method = "spearman") 
colnames(data.loess.sp) <- make.names(names, unique = TRUE)
rownames(data.loess.sp) <- make.names(names, unique = TRUE)

data.sp <- cor(data.norm, use = "complete.obs", method = "spearman") 
rownames(data.sp) <- make.names(names, unique = TRUE)
colnames(data.sp) <- make.names(names, unique = TRUE)

pairs(data.loess.sp, main = "Scatter plot Matrix - Loess Normalization")
pairs(data.sp, main = "Scatter plot Matrix - Global Median Normalization")


# Compare normalizations to quantile normalized data
array1.diff2 <- (maRf(dat[ , 1])) - (maRb(dat[ , 1]))
array2.diff2 <- (maRf(dat[ , 2])) - (maRb(dat[ , 2]))
array3.diff2 <- (maRf(dat[ , 3])) - (maRb(dat[ , 3]))
array4.diff2 <- (maRf(dat[ , 4])) - (maRb(dat[ , 4]))

samples <- cbind(array1.diff2, array2.diff2, array3.diff2, array4.diff2)
colnames(samples) <- c("Sample 1", "Sample 2", "Sample 3", "Sample 4")
samples.sort <- apply(samples, 2, sort) 
row.means <- rowMeans(samples.sort, na.rm = FALSE)

samp.order1 <- order(array1.diff2)  
samp.order2 <- order(array2.diff2)
samp.order3 <- order(array3.diff2)
samp.order4 <- order(array4.diff2)

samp.norm1 <- rep(NA, nrow(samples))
samp.norm2 <- rep(NA, nrow(samples))
samp.norm3 <- rep(NA, nrow(samples))
samp.norm4 <- rep(NA, nrow(samples))

samp.norm1[samp.order1] <- row.means
samp.norm2[samp.order2] <- row.means
samp.norm3[samp.order3] <- row.means
samp.norm4[samp.order4] <- row.means

quant.norm <- cbind(samp.norm1, samp.norm2, samp.norm3 samp.norm4)
head(samples)
head(quant.norm)


# Log (base 2) the new matrix and find Spearman's rank correlation
f.parse <- function(path=pa,file=fi,out=out.fi) {
	d <- read.table(paste(path,file,sep=""),skip=11,sep=",",header=T)
	u <- as.character(unique(d$Name))
	u <- u[u!=""]; u <- u[!is.na(u)];
	ref <- unique(as.character(d$Name[d$Type=="Reference"]))
	u <- unique(c(ref,u))
	hg <- c("B-actin","GAPDH","18S")
	hg <- toupper(hg)
	p <- unique(toupper(as.character(d$Name.1)))
	p <- sort(setdiff(p,c("",hg)))

	mat <- matrix(0,nrow=length(u),ncol=length(p))
	dimnames(mat) <- list(u,p)
	for (i in 1:length(u)) {
		print(paste(i,": ",u[i],sep=""))
		tmp <- d[d$Name %in% u[i],c(1:3,6,9)]
		g <- toupper(unique(as.character(tmp$Name.1)))
		g <- sort(setdiff(g,c("",hg)))
		
for (j in 1:length(g)) {
			v <- tmp[toupper(as.character(tmp$Name.1)) %in% g[j],5]
			v <- v[v!=999]
			v <- v[((v/mean(v))<1.5) & ((v/mean(v))>0.67)]	#gene j vector

			hv3 <- NULL
			for (k in 1:length(hg)) {	#housekeeping gene vector (each filtered by reps)
				hv <- tmp[toupper(as.character(tmp$Name.1)) %in% hg[k],5]
				hv <- hv[hv!=999]
				hv3 <- c(hv3,hv[((hv/mean(hv))<1.5) & ((hv/mean(hv))>0.67)]) 	
			}

			
			sv <- mean(as.numeric(v)) - mean(as.numeric(hv3))	#scaled value for gene j
			
			if(i==1) { #reference sample only
				mat[u[i],g[j]] <- sv
				next
			}
			
			mat[u[i],g[j]] <- sv - mat[u[1],g[j]]
		}
	}

	mat[1,][!is.na(mat[1,])] <- 0
	fc <- 2^(-1 * mat)
	write.table(t(c("Subject",dimnames(mat)[[2]])),paste(path,out,sep=""),quote=F,sep="\t",col.names=F,row.names=F)
	write.table(round(fc,3),paste(path,out,sep=""),quote=F,sep="\t",append=T,col.names=F)
}

pa <- "C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\HW2\\"
fi <- "Inflammation_qRT-PCR.csv"
out.fi <- "fold_chg_matrix.txt"

f.parse(pa,fi,out.fi)
