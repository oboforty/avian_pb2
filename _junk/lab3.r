install.packages("rJava")
install.packages("rmcfs")

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_131')
library(rmcfs)

###########################################################

alizadeh <- read.table("alizadeh.csv")

n = 10
alizadeh[1:n,c(1:n, length(alizadeh))]
alizadeh$class


n_projections = 3000


if(file.exists("result.rds")) {
  result <- readRDS("result.rds")
} else {
  
  result <- mcfs(class~., alizadeh, projections=n_projections,
                 projectionSize=0.1, splits=5, splitSetSize=0.66,
                 cutoffPermutations = 6, threadsNumber = 8)
  
  saveRDS(result, "result.rds")
}

# Values & Number of attributes:
head(result$RI)
dim(result$RI)

plot(result, type="distances")


# Task 5:
# check less significance values:
result$RI


# most significant ones:
result2 <- result$RI[1:result$cutoff_value,]
result2




# AFTER CUTOFF - Values & Number of attributes:
# for highest ranking genes:
head(result2)
dim(result2)




# Task 7:

gid <- build.idgraph(result, size = 20)
plot.idgraph(gid, label_dist = 0.3)
