library(devtools)
library(rmcfs)
library(dplyr)
library(R.ROSETTA)
library(VisuNet)


###############################
#                             #
#     Input Parameters        #
#                             #
###############################

filename <- "data/Project5.csv"
output_dir <- "output"

# MCFS parameters
n_projections <- 10000
proj_size <- 'auto'
cutoff_method <- 'kmeans'
cutoff_pe <- 20

splits <- 5
splitset_size <- 0.66

# Rosetta parameters
classifiers <- c("StandardVoter", "ObjectTrackingVoter", "NaiveBayesClassifier")
reducers <- c("Johnson", "Genetic")

cvNum <- 10
JohnsonParam <- list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, 
  Precompute=FALSE, Approximate=TRUE, Fraction=0.95
)
GeneticParam <- list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9,
  Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple")

# ROC parameters
host_clroc <- 'human'

###############################
#                             #
#        Load dataset         #
#                             #
###############################

load_protein_IS <- function(filename) {
  # loads file without header
  # this fixes a bug where T and F proteins are interpreted as booleans
  data <- read.table(filename, sep="\t", header=FALSE)

  # drop the header row and readd it as column names
  first_row <- data[1,]
  data <- data[-1,]
  colnames(data) <- as.character(unlist(first_row))
  data <- droplevels(data)

  return(data)
}

data = load_protein_IS(filename)


# Data metrics:
dim(data)
attributes(data)

attr = attributes(data)$names
#print(paste("Attributes per tree:", length(attr)*proj_size))

# Check if the dataset is balanced
table(data$Host)

###############################
#                             #
#            MCFS             #
#                             #
###############################

result <- mcfs(Host~., data, projections=n_projections,projectionSize=proj_size, splits=splits, splitSetSize=splitSetSize,
               cutoffPermutations = cutoff_pe, threadsNumber = 8)

head(result$RI)
plot(result, type="distances")

most_sig <- result$RI[1:result$cutoff_value,]
most_sig_names <- most_sig$attribute
rule_df <- select(data,most_sig_names,Host)

write.csv(rule_df,file = "sig_feat_table.csv")


#classifier = ObjectTrackingVoter. StandardVoter, NaiveBayesClassifier
#reducer = Johnson, Genetic


rules =rosetta(rule_df,roc = TRUE, clroc = "Human",classifier = "StandardVoter", reducer = "Johnson")
plotMeanROC(rules)

