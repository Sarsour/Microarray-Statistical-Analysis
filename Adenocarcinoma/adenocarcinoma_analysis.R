# Import data
dat <- read.table("lung_cancer.txt", header = T, row.names = 1)


# Load MASS Library and bind class names
library(MASS)
names <- colnames(dat)
dat <- data.frame(names, t(dat))
print(dim(dat))


# Create 2 different data sets
training <- rbind(dat[1:6, ], dat[11:16, ], dat[20:22, ])
test <- dat[ !(rownames(dat) %in% rownames(training)), ]
newtest.names <- test$names
test.names <- factor(gsub('[[:digit:]]+', '', newtest.names))
test$names <- NULL


# Run a classifier
newtraining.names <- training$names
train.names <- factor(gsub('[[:digit:]]+', '', newtraining.names))
training$names <- NULL
train.lda.2 <- lda(train.names ~ ., data = training[, c(1, 2)])
train.pred.2.out <- predict(train.lda.2, test[, c(1, 2)])
print(table(train.pred.2.out$class, test.names))


# Plot discriminant functions (first 2)
plot(range(train.pred.2.out$x[, 1]), range(train.pred.2.out$x[, 2]), type = "n", xlab = "LD1",ylab = "LD2", main = "LDA Plot - Training Set",)
points(train.pred.2.out$x, col = c(rep("Blue", 4), rep("Pink", 3), rep("Orange", 2)), pch = c(rep(16, 4), rep(17, 3), rep(18, 2)))
legend("bottomright", c("Adeno", "SCLC", "Normal"), col = c("Blue", "Pink", "Orange"), pch = c(16:18))


# Repeat with all
train.lda <- lda(train.names ~ ., data = training)
train.out <- predict(train.lda, test)
table(train.out$class, test.names)

plot(range(train.out$x[, 1]), range(train.out$x[, 2]), type = "n", xlab = "LD1", ylab = "LD2", main = "LDA Plot - Training Set - All Genes",)
points(train.out$x, col = c(rep("Blue", 4), rep("Pink", 3), rep("Orange", 2)), pch = c(rep(16, 4), rep(17, 3), rep(18, 2)))
legend("bottomright", c("Adeno", "SCLC", "Normal"), col = c("Blue", "Pink", "Orange"), pch = c(16:18))
