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
filename_results <- 'output/mcfs.rds'

# MCFS parameters
nr_projections <- 'auto'
projection_size <- 'auto'
cutoff_method <- 'permutations'
cutoff_permutations <- 20

splits <- 5
splitset_size <- 0.66

# Rosetta parameters
classifier <- "StandardVoter"
cvNum <- 10
reducer <- "Genetic"
JohnsonParam <- list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, 
  Precompute=FALSE, Approximate=TRUE, Fraction=0.95
)
GeneticParam <- list(
  Modulo=TRUE, BRT=FALSE, BRTprec=0.9,
  Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple")
underSample <- FALSE
underSampleNum <- 0
underSampleSize <- 0

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


dim(data)
attributes(data)
table(data$Host)


attr = attributes(data)$names
print(paste("Attributes per tree:", length(attr)*proj_size))



###############################
#                             #
#            MCFS             #
#                             #
###############################
?mcfs


###############################
#                             #
#          Rosetta            #
#                             #
###############################
?rosetta

ross_results <- rosetta(rule_df, discrete=TRUE, reducer=reducer, roc=TRUE, clroc=host_clroc, 
                        classifier=classifier, cvNum=cvNum, reducer=reducer, JohnsonParam=JohnsonParam, GeneticParam=GeneticParam, 
                        underSample=underSample, underSampleNum=underSampleNum, underSampleSize=underSampleSize)
