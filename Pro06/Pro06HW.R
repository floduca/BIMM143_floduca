#library(bio3d)
#s1 <- read.pdb("4AKE") # kinase with drug
#s2 <- read.pdb("1AKE") # kinase no drug
#s3 <- read.pdb("1E4Y") # kinase with drug
#s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
#s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
#s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
#s1.b <- s1.chainA$atom$b
#s2.b <- s2.chainA$atom$b
#s3.b <- s3.chainA$atom$b
#plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
#plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
#plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")


#Library calls the required library functions
library(bio3d)
#Protein plot takes an input x, where x is a pbd file
ProteinPlot <- function(x) {
  #Read pbd outputs a list when given a pdb file
  tablex <- read.pdb(x)
  #Trim PDB trims down the large pbd object into a smaller one only containing a subset of the original atoms
  chainx <- trim.pdb(tablex, chain = "A", elety = "CA")
  #Calls abd assigns a specific value in the trimmed list
  xb <- chainx$atom$b
  #Outputs a graph
  plotb3(xb, sse = chainx, typ = "l", ylab = "bfactor")
}
#Protein vector vectorizes all desired inputs to the Protein Plot function
ProteinVector = c("4AKE", "1AKE", "1E4Y")
#For loops through all the inputs in the vector
for (i in ProteinVector){
  ProteinPlot(i)
}