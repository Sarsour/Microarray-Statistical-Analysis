# Importing data to table
table <- read.table(file="C:\\Users\\Nidal\\Downloads\\sle_b_cell.txt", header=T, row.names=1)
print(colnames(table))

# Plotting SLE B cell samples vs Normal B cell samples
plot(table$sle.2, table$control.1, xlab = 'Normal', ylab = 'SLE', main = 'SLE B cell sample vs. Normal B cell sample – all probesets')
grid(col = "gray")

# Isolate and plot first 20
plot(table$sle.2[1:20], table$control.1[1:20], xlab = 'Normal', ylab = 'SLE', main = 'SLE B cell sample vs. Normal B cell sample – all probesets', pch=15, col = 'blue')

# Gene Profile Plot
plot(range(1:26), range(table["211881_x_at",]), type = "n",
xlab = "Sample", ylab = "Intensity", 
main = "IGLJ3 - Gene Profile Plot")
lines(c(1:26), as.numeric(table["211881_x_at",]))
grid(col = "grey") 

# Boxplot
f <- c(rep("SLE",17),rep("Control",9))
y <- as.numeric(table["211881_x_at",])
boxplot(y ~ f, xlab = "Sample", ylab = "Intensity", main = "IGLJ3 - Gene Profile Boxplot")