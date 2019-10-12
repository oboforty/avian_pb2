install.packages("rJava")
install.packages("rmcfs")

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_131')

library(rmcfs)
library(dplyr)
library(R.ROSETTA)
###########

# PREPARE DATASET
data <- read.table("data/Project5.csv", sep="\t", header=FALSE)
first_row <- data[1,]
data <- data[-1,]
colnames(data) <- as.character(unlist(first_row))
data <- droplevels(data)

###########
# DATA Poking
dim(data)
attributes(data)
table(data$Host)

###########

# MCFS feature selection
###?mcfs


n_projections <- 10000
cutoff_pe <- 20
proj_size <- 0.1
splits <- 5
splitSetSize <- 0.66

#if(file.exists("result.rds")) {
#  result <- readRDS("result.rds")
#} else {
  
result <- mcfs(Host~., data, projections=n_projections,projectionSize=proj_size, splits=splits, splitSetSize=splitSetSize,
               cutoffPermutations = cutoff_pe, threadsNumber = 8)
  
#  saveRDS(result, "result.rds")
#}

head(result$RI)
plot(result, type="distances")

most_sig <- result$RI[1:result$cutoff_value,]
rule_df <- select(data,sig_featname)


#Rosetta Rule Building
ross_results <- rosetta(rule_df, discrete = TRUE, reducer = "Genetic")
rule_table_info <- ross_results$main
viewRules(rule_table_info)


#Visunet Visualisation
visunet(rule_table_info)
