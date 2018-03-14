# Load in data

dat <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab6\\agingStudy11FCortexAffy.txt", header=T, row.names=1)
dat.ann <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab6\\agingStudy1FCortexAffyAnn.txt", header=T, row.names=1)


# Prepare vectors (age and gender)

dat.gender <- dat

dat.age.sort <- dat.ann[sort.list(dat.ann[, 2]), ]
age <- paste(dimnames(dat.age.sort)[[1]], dat.age.sort[, 2], dat.age.sort[, 1], sep = ".") 	
dat.age <- dat[, age]


# Run T test

g.g <- c(1394,  1474,  1917,  2099,  2367,  2428, 2625,  3168,  3181,  3641,  3832,  4526,
4731,  4863,  6062,  6356,  6684,  6787,  6900,  7223,  7244,  7299,  8086,  8652,
8959,  9073,  9145,  9389, 10219, 11238, 11669, 11674, 11793)

g.a <- c(25, 302,  1847,  2324,  246,  2757, 3222, 3675,  4429,  4430,  4912,  5640, 5835, 5856,  6803,  7229,  7833,  8133, 8579,  8822,  8994, 10101, 11433, 12039, 12353,
12404, 12442, 67, 88, 100)


t.test.all.genes <- function(x,s1,s2) {
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	t.out <- t.test(x1, x2, alternative = "two.sided", var.equal = T)
	out <- as.numeric(t.out$p.value)
	return(out)
}

rawp.gender <- apply(dat.gender[g.g, ], 1, t.test.all.genes, s1 = dat.ann[, 1] == "M", s2 = dat.ann[, 1] == "F") 
rawp.age <- apply(dat.age[g.a, ], 1, t.test.all.genes, s1 = age.sort[, 2] < 50, s2 = age.sort[, 2] >= 50) 

holm.gender <- mt.rawp2adjp(rawp.gender, proc = c("Holmes"))
holm.age <- mt.rawp2adjp(rawp.age, proc = c("Holmes"))
p.holm.gender <- p.adjust(rawp.gender, method = "holm") 
p.holm.age <- p.adjust(rawp.age, method = "holm") 


# Sort p-values and plot

rawp.gender.sort <- sort(rawp.gender)
p.holm.gender.sort <- sort(p.holm.gender) 
gen1 <- as.matrix(rawp.gender.sort)
gen2 <- as.matrix(p.holm.gender.sort)
gen <- cbind(gen1, gen2) 
colnames(gen) <- c("Raw P-value", "Adjusted P-Value") 
matplot(gen, type = "b", pch = 1, col = 1:2, main = "Gender P-values", ylab = "P-values") 
legend("topleft", legend = colnames(gen), pch = 1, col = 1:2) 

rawp.age.sort <- sort(rawp.age)
p.holm.age.sort <- sort(p.holm.age) 
age1 <- as.matrix(rawp.age.sort)
age2 <- as.matrix(p.holm.age.sort)
age <- cbind(age1, age2) 
colnames(age) <- c("Raw P-value", "Adjusted P-Value") 
matplot(age, type = "b", pch = 1, col = 1:2, main = "Age P-values", ylab = "P-values")
legend("topleft", legend = colnames(gen), pch = 1, col = 1:2)


# Bonferroni method

bon.gender <- p.adjust(rawp.gender, method = "bonferroni") 
bon.gender.sort <- sort(bon.gender) 
gen1.b <- as.matrix(rawp.gender.sort)
gen2.b <- as.matrix(bon.gender.sort)
gen.b  <- cbind(gen1.b, gen2.b) 
colnames(gen.b) <- c("Raw P-value", "Adjusted P-Value") 
matplot(gen.b, type = "b", pch = 1, col = 1:2, main = "Gender P-values", ylab = "P-values") 
legend("topleft", legend = colnames(gen.b), pch = 1, col = 1:2) 

bon.age <- p.adjust(rawp.age, method = "bonferroni") 
bon.age.sort <- sort(bon.age) 
age1.b <- as.matrix(rawp.age.sort)
age2.b <- as.matrix(bon.age.sort)
age.b  <- cbind(age1.b, age2.b) 
colnames(age.b) <- c("Raw P-value", "Adjusted P-Value") 
matplot(age.b, type = "b", pch = 1, col = 1:2, main = "Age P-values", ylab = "P-values") 
legend("topleft", legend = colnames(age.b), pch = 1, col = 1:2) 


# Read in data

tcga.dat <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab6\\tcga_brca_fpkm.txt", header=T, row.names=1)
tcga.dat.ann <- read.table(file="C:\\Users\\Nidal\\Documents\\JHU\\Gene Expression Data Analysis and Visualization\\Lab6\\tcga_brca_fpkm_sam.txt", header=T, row.names=1, na.strings=c("","NA"), sep="\t")
print(colnames(tcga.dat))


# Grep by gene 'GATA3'
tcga.dat.gata3 <- as.numeric(tcga.dat[grep("^GATA3_", tcga.dat$ROI, value=FALSE),])
	 

# Create binary vector for upper 25%
tcga.dat.gata3 <- sort(tcga.dat.gata3, decreasing = T)
tcga.dat.gata3[1:30] <- 1
tcga.dat.gata3[31:119] <- 0

group <- tcga.dat.gata3


f <- survfit(Surv(time, status),type="kaplan-meier",data=group)
plot(f,lwd=2,xlab='Weeks',ylab='S_hat(t)',main='Kaplan-Meier Plot')
summary(f)
survdiff(Surv(time, status) ~ x, data = aml)
