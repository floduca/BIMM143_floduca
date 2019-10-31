#' ---
#' title: "Data Visualization"
#' author: "Frank Loduca"
#' output: github_document
#' ---

# Class 05 Data Visualization
x <- rnorm(1000)
mean(x)
sd(x)
summary(x)
boxplot(x)
hist(x)
#rug(x)


weight <- read.table("bimm143_05_rstats/bimm143_05_rstats/weight_chart.txt")
