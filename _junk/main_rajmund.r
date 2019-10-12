
install.packages("rJava")
install.packages("rmcfs")
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_131')
library(rmcfs)


install.packages("devtools")
library(devtools)
install_github("komorowskilab/R.ROSETTA", force=TRUE)
library(R.ROSETTA)

###########

# PREPARE DATASET
data <- read.table("data/Project5.csv", sep="\t", header=FALSE)
first_row = data[1,]
data = data[-1,]
colnames(data) <- as.character(unlist(first_row))
data = droplevels(data)

###########
# DATA Poking
dim(data)
attr = attributes(data)$names
table(data$Host)

###########

# MCFS feature selection
?mcfs


# no. of subsets:
n_projections = 1500

# number of trees generated per subset
# no. of classifiers (trees) := splits * n_projections
splits = 5

# no./percent of loci in a subset:
proj_size = 0.2

# cutoff permutations
cutoff_pe = 20

splitSetSize = 0.66

# -----------------------------

print(paste("Attributes per tree:", length(attr)*proj_size))


# MCFS:
result <- mcfs(Host~., data, projections=n_projections,
                 projectionSize=proj_size, splits=splits, splitSetSize=splitSetSize,
                 cutoffPermutations = cutoff_pe, threadsNumber = 8)
# Save and load the results
saveRDS(result, "result.rds")

if(file.exists("result.rds")) {
  result <- readRDS("result.rds")
}

# Most significant
most_sig <- result$RI[1:result$cutoff_value,]
sig_cols = most_sig$attribute
dim(most_sig)


# Distances graph (projections convergence)
plot(result, type="distances")

# Interdependency graph
gid <- build.idgraph(result, size = result$cutoff_value)
plot.idgraph(gid, label_dist = 0.3)

# ---------------------------------
# Classifier

# Create information system
sig_data = data[c("id", sig_cols, "Host")]

# just to be sure
stopifnot( nrow(sig_data) == nrow(data) )

?rosetta
#out = rosetta(autcon)
out_johnson = rosetta(sig_data, discrete=TRUE, reducer="Johnson")

out_genetic = rosetta(sig_data, discrete=TRUE, reducer="Genetic")

rules_J = out_johnson$main[which(out_johnson$main$pValue < 0.05),]
rules_G = out_genetic$main[which(out_genetic$main$pValue < 0.05),]

viewRules(rules_G)

