setwd('D:/Nano-Clusters')

#install.packages('rgl')
library(rgl)     # for interactive 3d plot
#install.packages('plot3D')
library(plot3D)  # for 3d plot

PhyD <- read.table(file='Data/Phys_Descrip.txt', sep=',',dec='.', header=TRUE)
PhyD$Config <- as.factor(PhyD$Config)
for (j in c(1:4,6,7)){    
  PhyD[,j] <- as.integer(PhyD[,j])
}  

Cl_XYZ <- read.csv(file = 'Data/Cluster_XYZ.csv', 
                   dec='.', col.names=1:72, header=FALSE)
Nobs <- nrow(PhyD)


Var_names <- colnames(PhyD)

#Store the numerical descriptors in a matrix
PhyD_Mat <- matrix(unlist(PhyD[1:15] ), nrow=1839, ncol=15)
#summary(PhyD_Mat)







entry_nb <- 1752  # Index of the selected entry 

N_atom <- PhyD$N[entry_nb]
Mat_XYZ <- t(matrix(as.double(Cl_XYZ[entry_nb, 1:(3*N_atom)]), 3, N_atom))
colnames(Mat_XYZ) <- c('X','Y','Z')

noquote(sprintf('Entry number %i', (entry_nb)))

for (l in c(1:4,6,7)){
  print(sprintf('%s = %i', Var_names[l], PhyD_Mat[entry_nb,l]),quote=FALSE)
}
for (l in c(5,8:15)){
  print(sprintf('%s = %f', Var_names[l], PhyD_Mat[entry_nb,l]),quote=FALSE)
}
print(sprintf('%s = %s', Var_names[16], PhyD$Config[entry_nb]),quote=FALSE)

noquote('Atomic coordinates (rows are the coordinates of each atom):')
print(Mat_XYZ)


plot3d(Mat_XYZ[,1],Mat_XYZ[,2],Mat_XYZ[,3], 
       type='s', col = rainbow(N_atom), size=15, 
       xlab='X', ylab='Y', zlab='Z',
       ann = FALSE, axes = FALSE)
box3d()
rglwidget(elementId = "plot3drgl")





Cluster_FP_X1 <- function(Mat_XYZ){
  
  Barycenter <- apply(Mat_XYZ, 2, mean)
  Cent_Coo <- t(t(Mat_XYZ) - Barycenter)    
  
  CoMet = t(Cent_Coo) %*% Cent_Coo  
  Eig = eigen(CoMet) 
  
  return(c(Barycenter, sqrt(abs(Eig$values)), as.double(Eig$vectors)))
  
}  


entry_nb <- 1800
N_atom <- PhyD$N[entry_nb]
Mat_XYZ <- t(matrix(as.double(Cl_XYZ[entry_nb, 1:(3*N_atom)]), 3, N_atom))

noquote('Atomic coordinate:')
print(Mat_XYZ)
noquote('Cluster_FP_X1-fingerprint: ')
Cluster_FP_X1(Mat_XYZ)

ExtraFs <- matrix(nrow=Nobs, ncol=15)
for(i in 1:Nobs){
  N_atom <- PhyD$N[i]
  Mat_XYZ <- t(matrix(as.double(Cl_XYZ[i, 1:(3*N_atom)]), 3, N_atom))
  ExtraFs[i,] <- Cluster_FP_X1(Mat_XYZ)
}

write.csv(ExtraFs, 'StruFg.csv', row.names = F)
