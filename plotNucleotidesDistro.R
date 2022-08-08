args <- commandArgs(trailingOnly = TRUE)
fileName = args[1]
outFile = args[2]
outFileType = args[3]
plotTitle = args[4]

nData <- read.table(fileName, sep = "\t", header = TRUE) 

#import GGPLOT2
library(ggplot2)

#convert sequence length to factors and not integers
nData$length <- as.factor(nData$length)
#outFileName = paste(outFile, outFileType, sep=".")
#choice of output file type
if(outFileType == "png"){
        png(outFile)
}else if(outFileType == "pdf"){
        pdf(outFile)
}

nucPlot <- ggplot(nData, aes(fill=nucleotide, y=abundance, x=length)) + 
  geom_bar(position="fill", stat="identity") +
  xlab("Length (nt)") + 
  ylab("Frequency of first base") +
  ggtitle(plotTitle)

print(nucPlot)
dev.off()
