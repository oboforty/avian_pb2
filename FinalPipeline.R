library(devtools)
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_131')
library(rmcfs)
library(dplyr)
library(R.ROSETTA)
library(VisuNet)



##### FUNCTION DEFINITIONS #####

load_protein_IS <- function(filename) {
  
  #This function loads file without header and provides a work around for a bug where T and F proteins are interpreted as booleans
  data <- read.table(filename, sep="\t", header=FALSE)
  
  # drop the header row and readd it as column names
  first_row <- data[1,]
  data <- data[-1,]
  colnames(data) <- as.character(unlist(first_row))
  data <- droplevels(data)
  
  return(data)
}


load_mcfs_results <- function() {
  
  #This function loads saved mcfs result files from working directory
  if (exists("mcfs_result"))
    return (mcfs_result)
  
  if (file.exists(mcfs_result_filename)) {
    mcfs_result <- readRDS(mcfs_result_file)
    
    return (mcfs_result)
  }
}





##### FILE LOADING #####
filename <- "data/Project5.csv"
output_dir <- "output"
#mcfs_result_file <- "temp/mcfs_run.rds"   #Uncomment if using saved mcfs files
#features_filename <- "temp/most_sig_names.rds" #Uncomment if using saved features file


data = load_protein_IS(filename)

# Data metrics:
dim(data)
attributes(data)
attr = attributes(data)$names
#print(paste("Attributes per tree:", length(attr)*proj_size))

# Check if the dataset is balanced
table(data$Host)




##### MONTE CARLO FEATURE SELECTION #####

#Parameters
n_projections <- 10000
proj_size <- 'auto'
cutoff_method <- 'kmeans'
cutoff_pe <- 20
splits <- 5
splitset_size <- 0.66


mcfs_result <- mcfs(Host~., data, projections=n_projections,projectionSize=proj_size, splits=splits, splitSetSize=splitset_size,cutoffMethod = cutoff_method,
               cutoffPermutations = cutoff_pe, threadsNumber = 8)


head(mcfs_result$RI)
plot(mcfs_result, type="distances")


# Cache the results:
saveRDS(mcfs_result, mcfs_result_file)
mcfs_result <- load_mcfs_results()


# Get most significant names
most_sig <- mcfs_result$RI[1:mcfs_result$cutoff_value,]
most_sig_names <- most_sig$attribute


#saveRDS(most_sig_names, features_filename)
#most_sig_names <- readRDS(features_filename)


#Significant Feature Extraction
rule_df <- select(data, most_sig_names, Host)

#write.csv(rule_df,file = "sig_feat_table.csv")  #Uncomment to save the file as csv




##### ROSETTA RULE BASED CLASSIFICATION #####
classifier <- "StandardVoter" # StandardVoter, ObjectTrackingVoter, NaiveBayesClassifier
reducer <- "Johnson" # Johnson, Genetic

#JohnsonParam <- list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95)  #Uncomment if required
#GeneticParam <- list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple") #Uncomment if required

# Run Rosetta with ROC enabled
rules <- rosetta(rule_df,roc = TRUE, discrete=TRUE, clroc = "Human", classifier = classifier, reducer = reducer)
ruletable <- rules$main
viewRules(ruletable)

# Rule Quality
rules$quality

#Performance Visualisation
rules$ROC.stats
plotMeanROC(rules)



##### VISUNET VISUALISATION #####
visunet(ruletable)





