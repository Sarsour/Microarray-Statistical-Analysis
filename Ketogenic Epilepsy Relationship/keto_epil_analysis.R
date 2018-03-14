# Read data set
dat <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab5\\rat_KD.txt", header=T, row.names=1, na.strings="NA", blank.lines.skip=F)
print(dim(dat))


# Log2 the data and use T-test
dat.log2 <- log2(dat)
print(names(dat))
cl <- as.character(names(dat.log2))
dat.log2 <- dat.log2[ ,cl]

dat.control <- cl[1:6]
dat.keto <- cl[7:11]

t.test.all.genes <- function(x,s1,s2) {
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
	out <- as.numeric(t.out$p.value)
	return(out)
}

pv <- apply(dat.log2, 1, t.test.all.genes, s1 = dat.control, s2 = dat.keto)


# Histogram of p-values 
hist(pv,col="lightblue",xlab="p-values",main="P-value distance between\nControl and Ketogenic groups",cex.main=0.9)
abline(v=.05,col=2,lwd=2)

less.than.05 <- sum(pv < 0.05)
less.than.01 <- sum(pv < 0.01)
cons.p.value <- sum(pv < (0.05/length(pv)))

print(less.than.05)
print(less.than.01)
print(cons.p.value)


# Mean for each gene and fold change between groups
dat.control.m <- apply(dat[,dat.control],1,mean,na.rm=T)
dat.keto.m <- apply(dat[,dat.keto],1,mean,na.rm=T)
fold <- dat.control.m - dat.keto.m


# Min and max fold change value on linear scale and p-value threshold less than the Bonferroni threshold calculated
max.fold.change <- 2^max(fold)
min.fold.change <- 2^min(fold)
print(max.fold.change)
print(min.fold.change)

p.val.q6 <- names(pv[pv < cons.p.value & abs(2^fold) > 2])

print(p.val.q6)


# Transform p-value and volcano plot
p.trans <- -1 * log10(pv)

plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\nControl and Ketogenic group differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(2))],fold[(p.trans> -log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(2))],fold[(p.trans> -log10(.05)&fold< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))