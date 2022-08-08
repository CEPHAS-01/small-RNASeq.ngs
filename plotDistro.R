args <- commandArgs(trailingOnly = TRUE)
fileName = args[1]
minSize = args[2]
maxSize = args[3]
rows = args[4]
column = args[5]
outFile = args[6]
outFileType = args[7]

sizeDistro <- read.table(fileName, sep = "\t", header = T)

#import ggplot library
library(ggplot2)

#convert the sizes into factors and not integers
sizeDistro$Length = as.factor(sizeDistro$Length)

rows <- as.integer(rows)
column <- as.integer(column)

#plot the distribution
plotTitle <- paste("Length distribution of", minSize, "-", maxSize, "nt sequence across samples")


#outFileName = paste(outFile, outFileType, sep=".")
#choice of output file type
if(outFileType == "png"){
	png(outFile)
}else if(outFileType == "pdf"){
	pdf(outFile)
}

totalNumSamples <- rows * column

if(totalNumSamples > 6)
{
 distPlot <- ggplot(sizeDistro, aes( x = sizeDistro$Length, y = sizeDistro$Count, col = sizeDistro$Length)) +
  geom_bar( stat = "identity", fill="white") +
  xlab("Sequence length (nt)") +
  ylab("Count") + 
  ggtitle(plotTitle) +
  labs(col="Sequence Length") + 
  facet_wrap(~ sizeDistro$Sample, nrow = rows, ncol = column)  
}else{
 distPlot <- ggplot(sizeDistro, aes(x = sizeDistro$Length, y = sizeDistro$Count)) +
  geom_bar(
    aes(color = sizeDistro$Sample, fill = sizeDistro$Sample),
    stat = "identity", position = position_dodge(0.8),
    width = 0.5
  ) +
  xlab("Sequence length (nt)") +
  ylab("Count") +
  ggtitle(plotTitle) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}


print(distPlot)
dev.off()
