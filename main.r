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
data = droplevels(data)

###########
# DATA Poking
dim(data)
attributes(data)
table(data$Host)

###########

# MCFS feature selection
?mcfs


n_projections = 100
cutoff_pe = 20
proj_size = 0.1
splits = 5
splitSetSize = 0.66

#if(file.exists("result.rds")) {
#  result <- readRDS("result.rds")
#} else {
  
  result <- mcfs(Host~., data, projections=n_projections,
                 projectionSize=proj_size, splits=splits, splitSetSize=splitSetSize,
                 cutoffPermutations = cutoff_pe, threadsNumber = 8)
  
#  saveRDS(result, "result.rds")
#}

# Most significant
most_sig <- result$RI[1:result$cutoff_value,]


# Distances graph (projections convergence)
plot(result, type="distances")

# Interdependency graph
gid <- build.idgraph(result, size = result$cutoff_value)
plot.idgraph(gid, label_dist = 0.3)

