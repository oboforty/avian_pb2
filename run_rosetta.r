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

# Check if the dataset is balanced
table(data$Host)


# Reduce dataset with significant features (from MCFS runs)
#rule_df <- 



###############################
#                             #
#          Rosetta            #
#                             #
###############################

# running rosetta with changing parameters:
for (reducer in reducers) {
  for (classifier in classifiers) {
  
    ross_results <- rosetta(rule_df, reducer=reducer, classifier=classifier, clroc=host_clroc, 
                        roc=TRUE, discrete=TRUE, underSample=FALSE, cvNum=cvNum,
                        JohnsonParam=JohnsonParam, GeneticParam=GeneticParam)

    # change the results in the appropriate directory
    filename <- paste(output_dir, reducer, classifier, "result.rds", sep="/")
    result <- saveRDS(ross_results, filename)
  }
}
