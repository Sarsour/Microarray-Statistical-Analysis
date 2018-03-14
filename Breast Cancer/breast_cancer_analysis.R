# Load Data
dat <- read.table("sotiriou_data.txt", header = T, row.names = 1)
dat.ann <- read.table("Sotiriou_annotations.txt", header = T, row.names = 1)
print(dim(dat))


# PCA Plot
head(dat.ann[, "site"])
dat.pca <- prcomp(t(dat), cor=F)

plot(range(dat.pca$x[, 1]),range(dat.pca$x[, 2]), type = "n", xlab = "KIU", ylab = "OXF", main = "PCA Plot")
points(dat.pca$x[, 1][dat.ann$site == "KIU"], dat.pca$x[, 2][dat.ann$site == "KIU"], bg = "Red", pch = 21)
points(dat.pca$x[, 1][dat.ann$site == "OXF"], dat.pca$x[, 2][dat.ann$site == "OXF"], bg = "Blue", pch = 21)
legend("topright", c("KIU","OXF"), col = c("Red","Blue"), fill = c("Red", "Blue"))


# Scree Plot
dat.pca.var <- round(dat.pca$sdev^2 / sum(dat.pca$sdev^2) * 100, 2)
dat.pca.var

plot(c(1:length(dat.pca.var)), dat.pca.var, type = "b", xlab = "# Components", ylab = "% Variance", bg = "Blue", pch = 21)
title("Scree Plot")


# Classic and Nonmetric MDS PLots
	# Classic
dat.dist <- dist(t(dat))
dat.loc  <- cmdscale(dat.dist)
plot(dat.loc, type = "n")
points(dat.loc[, 1][dat.ann$site == "KIU"], dat.loc[, 2][dat.ann$site == "KIU"], col = "Red", pch=16, cex=1.5)
points(dat.loc[, 1][dat.ann$site == "OXF"], dat.loc[, 2][dat.ann$site == "OXF"], col = "Blue", pch=16, cex=1.5)
title(main="Classical MDS Plot")
legend("bottomright", c("KIU", "OXF"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
	# Nonmetric
library(MASS)
library(multtest)
dat.mds <- isoMDS(dat.dist)
plot(dat.mds$points, type = "n")
points(dat.mds$points[, 1][dat.ann$site == "KIU"], dat.mds$points[, 2][dat.ann$site == "KIU"], col = "Red", pch=16, cex=1.5)
points(dat.mds$points[, 1][dat.ann$site == "OXF"], dat.mds$points[, 2][dat.ann$site == "OXF"], col = "Blue", pch=16, cex=1.5)
title(main="Non-metric Plot")
legend("bottomright", c("KIU", "OXF"), col = c("Red", "Blue"), fill = c("Red", "Blue"))


# Weighted Graph Laplacian Plot
k.speClust2 <- function (X, qnt = NULL) {
	dist2full <- function(dis) {
		n <- attr(dis, "Size")
		full <- matrix(0, n, n)
		full[lower.tri(full)] <- dis
		full + t(full)
	}

	dat.dis <- dist(t(X),"euc")^2
	if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
	if(is.null(qnt)) {eps <- min(dat.dis[dat.dis!=0])}
	kernal <- exp(-1 * dat.dis/(eps))
	K1 <- dist2full(kernal)
	diag(K1) <- 0
	D = matrix(0,ncol=ncol(K1),nrow=ncol(K1))
	tmpe <- apply(K1,1,sum)
	tmpe[tmpe>0] <- 1/sqrt(tmpe[tmpe>0])
	tmpe[tmpe<0] <- 0
	diag(D) <- tmpe
	L <- D%*% K1 %*% D
	X <- svd(L)$u
	Y <- X / sqrt(apply(X^2,1,sum))
}

temp <- t(dat)
temp <- scale(temp, center = T, scale = T)
phi  <- k.speClust2(t(temp), qnt = NULL)

plot(
	range(phi[, 1]), range(phi[, 2]), 
	xlab = "Phi 1", ylab = "Phi 2", 
	main="Weighted Graph Laplacian Plot"
	)
points(phi[, 1][dat.ann$site == "KIU"], phi[, 2][dat.ann$site == "KIU"], col = "Red", pch = 16, cex = 1.5)
points(phi[, 1][dat.ann$site == "OXF"], phi[, 2][dat.ann$site == "OXF"], col = "Blue", pch = 16, cex = 1.5)
legend("bottomright", c("KIU", "OXF"), col = c("Red", "Blue"), fill = c("Red", "Blue"))
