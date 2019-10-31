Pro08 Handson
================
Paul Loduca
10/24/2019

``` r
fooddat <-  read.csv("UK_foods.csv", row.names = 1)
head(fooddat)
```

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

``` r
fooddat
```

    ##                     England Wales Scotland N.Ireland
    ## Cheese                  105   103      103        66
    ## Carcass_meat            245   227      242       267
    ## Other_meat              685   803      750       586
    ## Fish                    147   160      122        93
    ## Fats_and_oils           193   235      184       209
    ## Sugars                  156   175      147       139
    ## Fresh_potatoes          720   874      566      1033
    ## Fresh_Veg               253   265      171       143
    ## Other_Veg               488   570      418       355
    ## Processed_potatoes      198   203      220       187
    ## Processed_Veg           360   365      337       334
    ## Fresh_fruit            1102  1137      957       674
    ## Cereals                1472  1582     1462      1494
    ## Beverages                57    73       53        47
    ## Soft_drinks            1374  1256     1572      1506
    ## Alcoholic_drinks        375   475      458       135
    ## Confectionery            54    64       62        41

``` r
barplot(as.matrix(fooddat, beside = T, col = rainbow(nrow(fooddat))))
```

![](Pro08-Handson_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
pairs(fooddat, col=rainbow(10), pch=16)
```

![](Pro08-Handson_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
pca <-  prcomp(t(fooddat))
summary(pca)
```

    ## Importance of components:
    ##                             PC1      PC2      PC3       PC4
    ## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
    ## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    ## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

``` r
plot(pca$x[,1], pca$x[,2], xlab = "PC1", ylab = "PC2", xlim=c(-270,500))
text(pca$x[,1], pca$x[,2], colnames(fooddat), col = c("yellow", "red", "blue", "green"))
```

![](Pro08-Handson_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->
