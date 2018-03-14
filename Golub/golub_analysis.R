# Import libraries
library(Biobase)
library(multtest)
library(annotate)
data(golub)
library(limma)

# Matrix to data frame
dat <- as.data.frame(golub)
dimnames(dat)[[1]] <- paste("g", dimnames(dat)[[1]], sep = "")


# Set labels
dat.ann <- golub.cl
dimnames(dat)[[2]] <- ann 


# Wilcox test
wilcox.test.all.genes <- function(x,s1,s2) {
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	w.out <- wilcox.test(x1, x2, exact = FALSE, alternative = "two.sided", correct = TRUE, var.equal = T)
	out <- as.numeric(w.out$statistic)
	return(out)
}
original.wmw.run <- apply(dat, 1, wilcox.test.all.genes, s1 = colnames(dat) == 0, s2 = colnames(dat) == 1)


# WMV Test on multiple iterations
wilcox.iterate <- function(x){
	list <- list(1:500)
	for (i in 1:500){
		colnames(x) <- sample(colnames(x))
		stat <- apply(x, 1, wilcox.test.all.genes, s1 = colnames(x) == 0, s2 = colnames(x) == 1)
		list[i] <- max(stat)
	}
	result <- list
	return(result)
}

ptm <- proc.time()
stat.max <- wilcox.iterate(dat)
proc.time() - ptm
stat.max <- t(stat.max)
stat.max <- as.character(stat.max)


# Display 95% value test stat
cutoff <- sort(stat.max)[0.95 * length(stat.max)]
subset <- original.wmw.run[original.wmw.run > cutoff]
print(subset)


# Empirical Bayes method to compare
zero <- dat[colnames(dat) == 0]
one <- dat[colnames(dat) == 1]
x <- cbind(a = 1, b = c(rep(1, length(zero)), rep(0, length(one))))
fit <- lmFit(dat, x)
fit <- eBayes(fit)
topTable(fit)
attributes(fit)
pval <- fit$p.value[ , 2]


# Sort p-values and find commonalities
n <- length(subset)
sort <- sort(pval)
sort <- sort[1:n]
inter <- intersect(names(sort), names(subset))
length(inter)
print(inter)


#Compare with Student's T-test and graph p-values less than 0.01
t.run <- apply(dat, 1, t.test.all.genes, d1 = colnames(dat) == 0, d2 = colnames(dat) == 1)
t.run.p <- t.run[t.run < 0.01]
t.run.p <- as.matrix(t.run.p)
bayes.names <- as.matrix(fit$p.value)
bayes.names <- bayes.names[ , -2]
bayes.names <- as.matrix(bayes.names)
overlap <- merge(t.run.p, bayes.names, by = "row.names", all = FALSE)
overlap <- as.matrix(overlap)
rownames(overlap) <- overlap[ , 1]
overlap <- overlap[ , -1]
colnames(overlap) <- c("Student T-Test", "Empirical Bayes")
class(overlap) <- "numeric"

plot(c(1, nrow(overlap)), range(overlap), type = "n", xaxt ="n", ylab = "P-Value", xlab="Genes")
points(1:nrow(overlap), col = "Red", overlap[ , 2])
points(1:nrow(overlap), col = "Blue",  overlap[ , 1])
title(main="P-Value Plot of Student T-Test Vs. Empirical Bayes")
legend(300, 1 , colnames(overlap), col = c("Red", "Blue"), pch = 15, cex = 0.9)