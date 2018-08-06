### Data Generation for Scale Free Networks

# The following function counts shared edges between two graphs provided as 
# adjacency matrices
count_shared_edges <- function(adjTrue1, adjTrue2) {
  size = dim(adjTrue1)
  count_shared = 0
  for(i in 1:size[1])
  {
    for(j in 1:size[2]){
      if (adjTrue1[i,j]==adjTrue2[i,j] & adjTrue1[i,j]==1){
        count_shared = count_shared+1
      }
    }
  }
  return(count_shared)
}

library(igraph)
##############################################################
# The following code is repeated twice to acquire test and test2, 
# matrices which are adjusted to be "true" precision matrices for
# generated data

# Generate scale free network according to the Barabasi-Albert Model
graph.pa.1<-sample_pa(40)

# To plot graph if desired
plot.igraph(graph.pa.1,
            layout=layout.fruchterman.reingold,
            vertex.size=10,         # sets size of the vertex, default is 15
            vertex.label.cex=.5,    # size of the vertex label
            edge.arrow.size=.5        # sets size of the arrow at the end of the edge
)


# The function get.adjacency() converts a graph into an adjacency matrix.
pa.adjacency<-get.adjacency(graph.pa.1)
pa.adjacency

adj = as.matrix(pa.adjacency)
diag.40 = diag(40)
#To convert adjacency matrix from igraph package into a "starter" precision matrix

#test = adj + t(adj) +diag.40
test2= adj + t(adj) +diag.40

p=40

# to make test and test2 positive definite following 
# approach of Danaher et al (2014)
cc=3.4
dd=.95
#AA=test*cc
AA=test2*cc
AA = (AA+t(AA))
AA[AA>0] <- 1
AA = AA-diag(diag(AA))+diag(p)*dd
AA = AA/(as.matrix(rep(1,p)))%*%t(1.4*rowSums(abs(AA)))
AA = (AA+t(AA))/2
AA = AA-diag(diag(AA))+diag(p)*dd
#Omega1=AA
Omega2=AA
library(matrixcalc)
is.positive.definite(Omega1)
is.positive.definite(Omega2)
##############################################################
## To generate 3 "similar" networks
## NOTE: permutations may need to be adjusted to attain similarity levels like those in the paper

Omega1
Omega2
perm=1:40
perm2 = 1:40
sample(40,5) #3 26 1 35 8   reg order 1 3 8 26 35
sample(40,4)  #21, 39, 18, 12 reg order 12, 18, 21, 39
# Here we adjust the permutation vectors for edge location from the regular
# order to those attained in the sample
# This provides two randomly permuted edge orders for "similar" networks
perm[1]=3
perm[3]=26
perm[8]=1
perm[26]=35
perm[35]=8
perm2[12]=21
perm2[18]=39
perm2[21]=18
perm2[39]=12
Omega1.2 = Omega1[perm,perm]
Omega1.3 = Omega1[perm2,perm2]
is.positive.definite(Omega1.2)
is.positive.definite(Omega1.3)
adjTrue1=1*(Omega1!=0)
adjTrue2 = 1*(Omega1.2!=0)
adjTrue3 = 1*(Omega1.3!=0)

count_shared_edges(adjTrue1, adjTrue2)
#[1] 92
sum(colSums(adjTrue1))
#[1] 118
sum(colSums(adjTrue2))
#[1] 118
92/118
#[1] 0.779661
count_shared_edges(adjTrue1, adjTrue3)
#[1] 110
110/118
#[1] 0.9322034 
#In this case, networks 1 and 3 are "similar" 
# as at least 90% of edges are shared
count_shared_edges(adjTrue2, adjTrue3)
#[1] 84
84/118

Omega2
perm=1:40
perm2 = 1:40
sample(40,4) #29 32 15 11 reg order 11 15 29 32
sample(40,5) #31 2 15 16 26 reg order 2 15 16 26 31
perm[11]=29
perm[15]=32
perm[29]=15
perm[32]=11
perm2[2]=31
perm2[15]=2
perm2[16]=15
perm2[26]=16
perm2[31]=26

Omega2.2 = Omega2[perm,perm]
Omega2.3 = Omega2[perm2,perm2]
is.positive.definite(Omega2.2)
is.positive.definite(Omega2.3)
adj1 = 1*(Omega2!=0)
adj2 = 1*(Omega2.2!=0)
adj3 = 1*(Omega2.3!=0)
count_shared_edges(adj1, adj2) #110
sum(colSums(adj1)) #118
sum(colSums(adj2)) #118
sum(colSums(adj3)) #118
count_shared_edges(adj1, adj3) #90
count_shared_edges(adj2, adj3) #84

output_filename1 <- paste('./A1_g_p40_similar.csv', sep = "")
write.csv(Omega1, output_filename1)
output_filename2 <- paste('./A2_g_p40_similar.csv', sep = "")
write.csv(Omega1.3, output_filename2)
output_filename3 <- paste('./A3_g_p40_similar.csv', sep = "")
write.csv(Omega1.2, output_filename3)

output_filename1 <- paste('./A1_m_p40_similar.csv', sep = "")
write.csv(Omega2, output_filename1)
output_filename2 <- paste('./A2_m_p40_similar.csv', sep = "")
write.csv(Omega2.2, output_filename2)
output_filename3 <- paste('./A3_m_p40_similar.csv', sep = "")
write.csv(Omega2.3, output_filename3)

### Note: be sure to remove row and column names if reading into Matlab

## The following code generates three different scale-free networks for 
## both platforms, i.e. as in Setting 2 from the paper
graph1 <- sample_pa(40)
graph2 <- sample_pa(40)
graph3 <- sample_pa(40)
adj1 <- get.adjacency(graph1)
adj2 <- get.adjacency(graph2)
adj3 <- get.adjacency(graph3)

adj1 = as.matrix(adj1)
adj2 = as.matrix(adj2)
adj3 = as.matrix(adj3)

diag.40 = diag(40)

#test = adj1 + t(adj1) + diag.40
#test = adj2 + t(adj2) + diag.40
test = adj3 + t(adj3) + diag.40
p=40
#to make positive definite
cc=3.4
dd=.95
AA=test*cc
AA = (AA+t(AA))
AA[AA>0] <- 1
AA = AA-diag(diag(AA))+diag(p)*dd
AA = AA/(as.matrix(rep(1,p)))%*%t(1.4*rowSums(abs(AA)))
AA = (AA+t(AA))/2
AA = AA-diag(diag(AA))+diag(p)*dd
#Omega1=AA
#Omega2 = AA
Omega3 = AA
is.positive.definite(Omega1)
is.positive.definite(Omega2)
is.positive.definite(Omega3)
adj1 = 1*(Omega1!=0)
adj2 = 1*(Omega2!=0)
adj3 = 1*(Omega3!=0)
count_shared_edges(adj1, adj2) #46
sum(colSums(adj1)) #118
sum(colSums(adj2)) #118
sum(colSums(adj3)) #118
count_shared_edges(adj1, adj3) #48
count_shared_edges(adj2, adj3) #44

output_filename1 <- paste('A1_g_p40_different.csv', sep = "")
write.csv(Omega1, output_filename1)
output_filename2 <- paste('A2_g_p40_different.csv', sep = "")
write.csv(Omega2, output_filename2)
output_filename3 <- paste('A3_g_p40_different.csv', sep = "")
write.csv(Omega3, output_filename3)

graph1 <- sample_pa(40)
graph2 <- sample_pa(40)
graph3 <- sample_pa(40)
adj1 <- get.adjacency(graph1)
adj2 <- get.adjacency(graph2)
adj3 <- get.adjacency(graph3)

adj1 = as.matrix(adj1)
adj2 = as.matrix(adj2)
adj3 = as.matrix(adj3)

diag.40 = diag(40)

test = adj1 + t(adj1) + diag.40
#test = adj2 + t(adj2) + diag.40
#test = adj3 + t(adj3) + diag.40
p=40
#to make positive definite
cc=3.4
dd=.95
AA=test*cc
AA = (AA+t(AA))
AA[AA>0] <- 1
AA = AA-diag(diag(AA))+diag(p)*dd
AA = AA/(as.matrix(rep(1,p)))%*%t(1.4*rowSums(abs(AA)))
AA = (AA+t(AA))/2
AA = AA-diag(diag(AA))+diag(p)*dd
Omega1=AA
#Omega2 = AA
#Omega3 = AA
is.positive.definite(Omega1)
is.positive.definite(Omega2)
is.positive.definite(Omega3)
adj1 = 1*(Omega1!=0)
adj2 = 1*(Omega2!=0)
adj3 = 1*(Omega3!=0)
count_shared_edges(adj1, adj2) #58
sum(colSums(adj1)) #118
sum(colSums(adj2)) #118
sum(colSums(adj3)) #118
count_shared_edges(adj1, adj3) #48
count_shared_edges(adj2, adj3) #48

output_filename1 <- paste('A1_m_p40_different.csv', sep = "")
write.csv(Omega1, output_filename1)
output_filename2 <- paste('A2_m_p40_different.csv', sep = "")
write.csv(Omega2, output_filename2)
output_filename3 <- paste('A3_m_p40_different.csv', sep = "")
write.csv(Omega3, output_filename3)