args <- commandArgs(trailingOnly = TRUE)
fileName = args[1]
outFile = args[2]
logFC = as.numeric(args[3])
pVal = as.numeric(args[4])
outFileType = args[5]

#volcano plot of expression data
vData <- read.table(fileName, sep = "\t", header = T)
#head(vData)

#load the library
library(EnhancedVolcano)

outFileName = paste(outFile,outFileType,sep=".")

if(outFileType == "png"){
	#outFileName = paste(outFile,outFileType,sep=".")
        png(outFileName)
}else if(outFileType == "pdf"){
        pdf(outFileName)
}

outPlot <- EnhancedVolcano(vData,
                lab = rownames(vData),
                x = 'logFC',
                y = 'pVal',
                pCutoff = pVal,
                FCcutoff = logFC,
                pointSize = 4.0,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0)


print(outPlot)
dev.off()
