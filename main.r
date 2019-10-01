# INSTALL


###########

# PREPARE DATASET
data <- read.table("data/Project5.csv", sep="\t", header=FALSE)
first_row = data[1,]
data = data[-1,]
colnames(data) <- as.character(unlist(first_row))


###########

# MCFS feature selection