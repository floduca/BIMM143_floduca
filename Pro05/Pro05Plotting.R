weight <- read.table("bimm143_05_rstats/bimm143_05_rstats/weight_chart.txt",
                     header = TRUE,
                     sep = "")
plot(weight,
     #Type sets the type of graph
     type = "o",
     #pch sets the symbol used
     pch = 15, 
     #cex is magnification
     cex = 1.5,
     #lwd is line width
     lwd = 2, 
     #ylim/xlim sets axis parameters
     ylim = c(2, 10), 
     #xlab/ylab sets labels of axis
     xlab = "Age (Months)",
     ylab = "Weight (Kg)",
     #Main sets the overall title
     main = "Weight")

mouse <- read.table("bimm143_05_rstats/bimm143_05_rstats/feature_counts.txt",
           header = TRUE,
           sep = "\t")
#mar sets the number of lines of a given margin mar = c(bottom, left, top, right)
par(mar = c(3.1, 11.1, 4.1, 2))
barplot(
  #Calls only the count column of Mouse
  height = mouse$Count,
  #Sets bars to horizontal
        horiz = TRUE,
  #Calls the feature column to name the bars
  names.arg = mouse$Feature,
  #Sets the main title
  main = "Mouse Genomic Features",
  #Sets the axis label to horizontal
  las = 1,)
#reading a thing
malefemale <- read.table("bimm143_05_rstats/bimm143_05_rstats/male_female_counts.txt",
                         header = TRUE,
                         sep = "\t")
#bargraph
barplot(height = malefemale$Count,
        col = c("blue2", "red2"))
#reading a thing
updown <- read.delim("bimm143_05_rstats/bimm143_05_rstats/up_down_expression.txt",
                     header = TRUE,
                     sep = "\t")
nrow(updown)
#Setting the color palette
palette(c("red", "black", "green"))
#Gets information about the State column of updown
table(updown$State)
plot(updown$Condition1, updown$Condition2,
     #Color by conditional state
     col = updown$State)
levels(updown$State)