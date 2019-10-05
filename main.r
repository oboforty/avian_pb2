#### Packages Required ###

# UNCOMMENT THESE LINES IF THE PACKAGES ARE NOT INSTALLED

# install.packages("rJava")
# install.packages("rmcfs")
# install.packages("dplyr")
# devtools::install_github("komorowskilab/VisuNet")

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_131')

library(rmcfs)
library(dplyr)
library(R.ROSETTA)
library(VisuNet)
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

n_projections <- 100
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


# Extract Features
most_sig <- result$RI[1:result$cutoff_value,]
sig_featname <- most_sig$attribute
rule_df <- select(data,sig_featname,Host)


#Interdependency Graph
gid <- build.idgraph(result, size = 20)
plot.idgraph(gid, label_dist = 0.3)


#Rosetta Rule Building
ross_results <- rosetta(rule_df, discrete = TRUE, reducer = "Genetic")
rule_table_info <- ross_results$main
viewRules(rule_table_info)


#ROC Plot
roc_rules <- rosetta(rule_df, discrete = TRUE, reducer = "Genetic", roc = TRUE, clroc = 'Avian')
plotMeanROC(roc_rules)


#Visunet Visualisation
visunet(rule_table_info)
