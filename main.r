install.packages("rJava")
install.packages("rmcfs")

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_131')
library(rmcfs)


###########

# PREPARE DATASET
data <- read.table("data/Project5.csv", sep="\t", header=FALSE)
first_row = data[1,]
data = data[-1,]
colnames(data) <- as.character(unlist(first_row))

###########
# DATA Poking
dim(data)
attributes(data)
table(data$Host)

###########

# MCFS feature selection
?mcfs


n_projections = 3000


if(file.exists("result.rds")) {
  result <- readRDS("result.rds")
} else {
  
  #result <- mcfs(...)
  
  saveRDS(result, "result.rds")
}
