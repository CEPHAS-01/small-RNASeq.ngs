args <- commandArgs(trailingOnly = TRUE)
file_in = args[1]   ### es."/biodata/ArabidopsisVirus/smallRNA/t-test/G3_I3_M3_N3_raw_mirna.txt"  input raw freq
control = args[2]  ###es "mock" 
exp = args[3] ####es."CaMV"  
output_name = args[4] ###es "mock_CaMV"  output file
n1 = as.integer(args[5])  ##es. 2  number of replicates of control
n2 = as.integer(args[6])  ##es. 2  number of replicates of experiment
pval = as.double(args[7]) ### es. 0.05  ## pvalue threshold
countpm = as.integer(args[8]) ### es. 1 ##keep only those genes that have at least countpm read per million in at least nsamples samples
nsamples = as.integer(args[9]) ### es. 2  ##keep only those genes that have at least countpm read per million in at least nsamples samples
#dval = args[8] ## es 25 [50/(n1+n2-2)]

#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("MASS")
library(MASS)
library(edgeR)
raw.data <- read.table( file = file_in , header = TRUE )
counts <- raw.data[ , -c(1,ncol(raw.data)+1) ]
#counts <- raw.data[ , -c(1,ncol(raw.data)) ]
rownames( counts ) <- raw.data[ , 1 ]
colnames<-paste(raw.data[1,],c(1:n1,1:n2),sep="")
#colnames( counts ) <- paste(c(rep(control,n1),rep(exp,n2)),c(1:n1,1:n2),sep="")
#head(counts)
#dim( counts )
#colSums( counts )         ## Library Sizes
#colSums( counts ) / 1e06     ##Library Sizes in millions of reads
#table( rowSums( counts ) )[ 1:30 ]  # Number of genes with low counts
group <- c(rep(control, n1) , rep(exp, n2))
cds <- DGEList( counts , group = group )
#names( cds )
#head(cds$counts)  # original count matrix
#cds$samples  # contains a summary of your samples
#sum( cds$all.zeros ) # How many genes have 0 counts across all samples
#cds # or type the name of the object
# keep only those genes that have at least 1 read per million in at least 2 samples
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > countpm ) >= nsamples , ]
#dim( cds )
cds <- calcNormFactors( cds )
#cds$samples
# effective library sizes
#cds$samples$lib.size * cds$samples$norm.factors
#title<-paste(control,exp, sep = "_")
bcv_title<-paste(output_name,"plotBCV.pdf",sep="_")
mds_title<-paste(output_name,"plotMDS.pdf",sep="_")
out_title<-paste(output_name,"DE.txt",sep="_")
#pdf( "plotMDS.pdf" , width = 7 , height = 7 )
pdf( mds_title , width = 7 , height = 7 )
# To view the plot immediately
plotMDS(cds, method="bcv",  main = "MDS Plot for Count Data", col=as.numeric(cds$samples$group))
legend("bottomleft", as.character(unique(cds$samples$group)), col=1:3, pch=20)

#plotMDS( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )
#crea pdf

dev.off()
# Estimating Dispersions
cds <- estimateCommonDisp( cds )
#names( cds )
# The estimate
#cds$common.dispersion

#cds <- estimateTagwiseDisp( cds , prior.df= 10 )
cds <- estimateTagwiseDisp( cds  )
#names( cds )
#summary( cds$tagwise.dispersion )
pdf( bcv_title , width = 7 , height = 7 )
#pdf( "plotBCV.pdf")
plotBCV(cds)
dev.off()
#meanVarPlot <- plotMeanVar( cds , show.raw.vars=TRUE ,
#show.tagwise.vars=TRUE ,
#show.binned.common.disp.vars=FALSE ,
#show.ave.raw.vars=FALSE ,
#NBline = TRUE ,
#nbins = 100 ,
#pch = 16 ,
#xlab ="Mean Expression (Log10 Scale)" ,
#ylab = "Variance (Log10 Scale)" ,
#main = "Mean-Variance Plot" )

de.tgw <- exactTest( cds , pair = c( control , exp ) )
#de.poi <- exactTest( cds , dispersion = 1e-06 , pair = c( "mock" , "CaMV" ) )
# Top tags for tagwise analysis
options( digits = 3 ) # print only 3 digits
#topTags( de.tgw , n = 20 , sort.by = "p.value" ) # top 20 DE genes
# Back to count matrix for tagwise analysis
#cds$counts[ rownames( topTags( de.tgw , n = 15 )$table ) , ]
# Sort tagwise results by Fold-Change instead of p-value
#resultsByFC.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) , sort.by = "logFC" )$table
#head( resultsByFC.tgw )

# Store full topTags results table
resultsTbl.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) )$table
#resultsTbl.poi <- topTags( de.poi , n = nrow( de.poi$table ) )$table
#head( resultsTbl.tgw )
#write.table( resultsTbl.tgw , file = file_out , sep="\t", row.names = TRUE )



colnames( resultsTbl.tgw ) <- c(  "logFC" , "logConc" , "pVal.Tgw" , "adj.pVal.Tgw" )
wh.rows.tgw <- match( rownames( resultsTbl.tgw ) , rownames( cds$counts ) )

# Change column names to be specific to the analysis, logConc and logFC are the same in both.
# Below provides the info to re-order the count matrix to be in line with the order of the results.

#head( wh.rows.tgw )
# Tagwise Results


combResults.tgw <- cbind( resultsTbl.tgw ,"Tgw.Disp" = cds$tagwise.dispersion[ wh.rows.tgw ] ,"UpDown.Tgw" = decideTestsDGE( de.tgw , p.value = pval )[ wh.rows.tgw ] ,cds$counts[ wh.rows.tgw , ] )
#head( combResults.tgw )
# Common Results
options( digits = 3 )
#write.matrix( combResults.tgw , file = out_title, sep = "\t" , row.names = TRUE)
combResults.tgw<-round(combResults.tgw, digits=6)
combResults.tgw<- cbind(ID = rownames(combResults.tgw), combResults.tgw)
rownames(combResults.tgw) <- NULL
#names(dimnames(combResults.tgw)) <- c("id", "")
write.table( combResults.tgw , file = out_title, sep = "\t"  , row.names=FALSE, col.names=TRUE, quote=FALSE )


