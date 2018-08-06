### Simulated Example of Metabolite Selection
library(corrplot)
library(grDevices)
## Generating correlated normal data 
## with mean = 0 and variance = 1 for simplicity

# Correlation matrix attained from rounding to three decimals
# the correlation of the 28th - 42nd metabolite abundances 
# matching to Kegg ID "C01194"
# of all subjects from Reg Auto pathway

corr.matrix = matrix(cbind(1.000, -0.239, -0.216, -0.169,  0.196,  0.201, -0.194, -0.163, -0.073, -0.092, -0.160, -0.001,  0.015,  0.285,  0.197,
                           -0.239,  1.000,  0.985,  0.213,  0.433,  0.453,  0.136,  0.122, -0.377,  0.569, -0.043, -0.012,  0.031,  0.147,  0.082,
                           -0.216,  0.985,  1.000,  0.217,  0.420,  0.434,  0.134,  0.135, -0.378,  0.560, -0.052, -0.005,  0.039,  0.164,  0.093,
                           -0.169,  0.213,  0.217,  1.000, -0.092, -0.098,  0.967,  0.952,  0.095, -0.096,  0.445,  0.262,  0.286, -0.163, -0.256,
                           0.196,  0.433,  0.420, -0.092,  1.000,  0.956, -0.149, -0.158, -0.410,  0.268, -0.161, -0.042, -0.009,  0.181,  0.197,
                           0.201,  0.453,  0.434, -0.098,  0.956,  1.000, -0.147, -0.181, -0.380,  0.274, -0.169, -0.022,  0.012,  0.194,  0.216,
                           -0.194,  0.136,  0.134,  0.967, -0.149, -0.147,  1.000,  0.935,  0.149, -0.141,  0.502,  0.276,  0.294, -0.229, -0.314,
                           -0.163,  0.122,  0.135,  0.952, -0.158, -0.181,  0.935,  1.000,  0.180, -0.143,  0.454,  0.281,  0.302, -0.192, -0.272,
                           -0.073, -0.377, -0.378,  0.095, -0.410, -0.380,  0.149,  0.180,  1.000, -0.274,  0.209,  0.158,  0.114, -0.164, -0.160,
                           -0.092,  0.569,  0.560, -0.096,  0.268,  0.274, -0.141, -0.143, -0.274,  1.000, -0.264, -0.204, -0.128,  0.193,  0.204,
                           -0.160, -0.043, -0.052,  0.445, -0.161, -0.169,  0.502,  0.454,  0.209, -0.264,  1.000,  0.350,  0.350, -0.217, -0.175,
                           -0.001, -0.012, -0.005,  0.262, -0.042, -0.022,  0.276,  0.281,  0.158, -0.204,  0.350,  1.000,  0.947, -0.016, -0.041,
                           0.015,  0.031,  0.039,  0.286, -0.009,  0.012,  0.294,  0.302,  0.114, -0.128,  0.350,  0.947,  1.000,  0.003, -0.028,
                           0.285,  0.147,  0.164, -0.163,  0.181,  0.194, -0.229, -0.192, -0.164,  0.193, -0.217, -0.016,  0.003,  1.000,  0.768,
                           0.197,  0.082,  0.093, -0.256,  0.197,  0.216, -0.314, -0.272, -0.160,  0.204, -0.175, -0.041, -0.028,  0.768,  1.000),nrow=15)

C = t(chol(corr.matrix))
nvars = dim(C)[1]
numobs = 20
set.seed(1)
rand.norm.vars = matrix(rnorm(nvars*numobs,0,1), nrow=nvars, ncol=numobs);
X = C %*% rand.norm.vars
newX = t(X)
cor.Data = as.data.frame(newX)
orig.Data = as.data.frame(t(rand.norm.vars))
#Names of "Metabolite Variables"
colnames(cor.Data) = c("i","ii","iii", "iv", "v",
                       "vi", "vii", "viii", "ix", "x",
                       "xi","xii","xiii","xiv","xv")
# Names of "Observations"
rownames(cor.Data) = c("A","B","C","D","E","F","G","H","I","J",
                       "K","L","M","N","O","P","Q","R","S","T")
cor(cor.Data)
plot(cor.Data)

plot(orig.Data)
generatedDataForExample <- t(cor.Data)
write.csv(generatedDataForExample, file = "corData_MetabSelectionExample.csv")

# Function to remove extraneous metabolite variants by extracting less correlated variants 
#until the first principle component analysis variant explains 98 % of variance 
# Arguments:
#     sub = list of subject measurements for duplicate kegg_id variants
# Return:
#     list of (removed, selected)
#     removed(less correlated subset), 
#     selected(remaining highly correlated subset where first principle component explains most of variance)

func1 <- function(sub){
  #Center and standardize input
  standardized.sub <- as.data.frame(scale(sub))
  # perform principal component analysis
  pca.sub <- prcomp(standardized.sub)
  # compute percent of variance contributed by each component
  perc.var.all <- (pca.sub$sdev)^2/sum(pca.sub$sdev^2)
  #while first component explains less than 98 %, 
  #select and remove into a subset the least correlated variable
  while (perc.var.all[1]<.98){
    new.removed <- as.data.frame(rbind(new.removed,sub[(which.min(abs(cor(t(sub))[1,]))),]))
    sub <- sub[-(which.min(abs(cor(t(sub))[1,]))),]
    standardized.subset <- as.data.frame(scale(sub))
    pca <- prcomp(standardized.subset)
    perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
  }
  # remaining more correlated variables where 1st pca variant 
  # contributes to at least 98 % of variance
  selected <- sub
  # set output to be less correlated subset and correlated subset
  output <- list(new.removed, selected)
  return(output)
}

## Note this function reads in the data with variables
## as rows and observations as columns
new.removed <- NULL
selected <- matrix()
# Input Data read in so that rownames output the same for all iterations of function
# otherwise some iterations give row index while first iteration gives metabolite names
input <- read.csv("corData_MetabSelectionExample.csv",header=TRUE)
#remove column of "metabolite names"
noRowName <- input[,2:21]
a <- func1(noRowName)
removed <- a[[1]]
rownames(removed)
# [1] "6"  "5"  "10" "8"  "13" "4"  "12" "11" "9"  "7"  "3"  "2"  "15"
selected <- a[[2]]
rownames(selected)
# [1] "1"  "14"    
# The above rows are remaining as a highly correlated subset where the variable
# correspoonding to row "1" explains more than 98 % of variance of subset
std.selected <- as.data.frame(scale(selected))
# The following commands illustrate percent of variance explained by pca of the 
# correlated extracted subset
pca <- prcomp(std.selected)
perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
perc.var.all
#[1] 1.000000e+00 3.882675e-32

a <- func1(removed)
removed2 <- a[[1]]
rownames(removed2)
# [1] "4"  "10" "8"  "7"  "12" "13" "11" "15" "9"  "2"  "3" 
selected2 <- a[[2]]
rownames(selected2)
# [1] "6" "5"
std.selected <- as.data.frame(scale(selected2))
pca <- prcomp(std.selected)
perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
perc.var.all
# [1] 1.000000e+00 1.232595e-32

a <- func1(removed2)
removed3 <- a[[1]]
rownames(removed3)
# [1] "9"  "3"  "2"  "10" "15" "12" "13" "11" "7" 
selected3 <- a[[2]]
rownames(selected3)
# [1] "4" "8"
std.selected <- as.data.frame(scale(selected3))
pca <- prcomp(std.selected)
perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
perc.var.all
# [1] 1.000000e+00 3.420452e-32

a <- func1(removed3)
removed4 <- a[[1]]
rownames(removed4)
# [1] "13" "7"  "12" "11" "15" "10" 
selected4 <- a[[2]]
rownames(selected4)
# [1] "9" "3" "2"
std.selected <- as.data.frame(scale(selected4))
pca <- prcomp(std.selected)
perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
perc.var.all
# [1] 9.851678e-01 1.483215e-02 1.567633e-33

a <- func1(removed4)
removed5 <- a[[1]]
rownames(removed5)
# [1] "15" "7"  "10" "11"
selected5 <- a[[2]]
rownames(selected5)
# [1] "13" "12"
std.selected <- as.data.frame(scale(selected5))
pca <- prcomp(std.selected)
perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
perc.var.all
# [1] 1.000000e+00 3.420452e-32

a <- func1(removed5)
removed6 <- a[[1]]
rownames(removed6)
# [1] 10" "11"
selected6 <- a[[2]]
rownames(selected6)
# [1] "15" "7"
std.selected <- as.data.frame(scale(selected6))
pca <- prcomp(std.selected)
perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
perc.var.all
# [1] 1.000000e+00 3.882675e-32

a <- func1(removed6)
removed7 <- a[[1]]
rownames(removed7)
# NULL
selected7 <- a[[2]]
rownames(selected7)
# [1] "10" "11"
std.selected <- as.data.frame(scale(selected7))
pca <- prcomp(std.selected)
perc.var.all <- (pca$sdev)^2/sum(pca$sdev^2)
perc.var.all
# [1] 1.000000e+00 1.32504e-32

## Selected Subset Rows are then the first of each "extracted" group
## 1,6,4,9,13,15,10

## Selected metabolites for network analysis would then be this subset of 
## metabolites and observations
input_subset <- input[c(1,6,4,9,13,15,10),2:21]
newInput = as.data.frame(t(input_subset))
colnames(newInput) = c("i","vi","iv","ix","xiii",
                       "xv", "x")

corNewSubset = cor(newInput)
plot(newInput)
