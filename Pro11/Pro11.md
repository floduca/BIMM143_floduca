Pro11
================
Paul Loduca
11/5/2019

\#\#PDB databsed for biomolecular structure data

``` r
PDBdata <- read.csv("Data Export Summary.csv")
PDBdata
```

    ##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other
    ## 1               X-Ray   131278          2059               6759     8
    ## 2                 NMR    11235          1303                261     8
    ## 3 Electron Microscopy     2899            32                999     0
    ## 4               Other      280             4                  6    13
    ## 5        Multi Method      144             5                  2     1
    ##    Total
    ## 1 140104
    ## 2  12807
    ## 3   3930
    ## 4    303
    ## 5    152

``` r
methodproportions <- PDBdata$Total / sum(PDBdata$Total)
methodproportions
```

    ## [1] 0.890702879 0.081419744 0.024984742 0.001926305 0.000966331

``` r
proportionprotein <- sum(PDBdata$Proteins)/sum(PDBdata$Total)
proportionprotein
```

    ## [1] 0.9271437

\#\#1HSG Structure

``` r
library(bio3d)
pdb <-  read.pdb("1hsg.pdb")
whole <- atom.select(pdb)
whole
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb)
    ## 
    ##    Atom Indices#: 1686  ($atom)
    ##    XYZ  Indices#: 5058  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
a.inds <- atom.select(pdb, chain = "A")
a.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, chain = "A")
    ## 
    ##    Atom Indices#: 801  ($atom)
    ##    XYZ  Indices#: 2403  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
ca.inds <-  atom.select(pdb, "calpha", chain = "A")
ca.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "calpha", chain = "A")
    ## 
    ##    Atom Indices#: 99  ($atom)
    ##    XYZ  Indices#: 297  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
cab.inds <- atom.select(pdb, elety =c("CA", "CB"), chain = "A", resno = 10:20)
cab.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, elety = c("CA", "CB"), chain = "A", 
    ##     resno = 10:20)
    ## 
    ##    Atom Indices#: 20  ($atom)
    ##    XYZ  Indices#: 60  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
protein.inds <- atom.select(pdb, string = "protein", value = TRUE)
write.pdb(protein.inds, file = "1hsg_protein")
ligand.inds <- atom.select(pdb, string = "ligand", value = TRUE)
write.pdb(ligand.inds, file = "1hsg_ligand")
ligand.inds
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call
