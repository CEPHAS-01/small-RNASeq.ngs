args <- commandArgs(trailingOnly = TRUE)
file_in = args[1]   ### es."/biodata/ArabidopsisVirus/smallRNA/t-test/G3_I3_M3_N3_raw_mirna.txt"  input raw freq
control = args[2]  ###es "mock" 
exp = args[3] ####es."CaMV"  
file_out = args[4] ###es "/biodata/ArabidopsisVirus/smallRNA/t-test/G3I3_M3N3_edger_mirna.txt"  output file
n1 = args[5]  ##es. 2  number of replicates of control
n2 = args[6]  ##es. 2  number of replicates of experiment
pval = args[7] ### es. 0.05  ## pvalue threshold
#dval = args[8] ## es 25 [50/(n1+n2-2)]

#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
raw.data <- read.table( file = file_in , header = TRUE )
counts <- raw.data[ , -c(1,ncol(raw.data)+1) ]
rownames( counts ) <- raw.data[ , 1 ]
colnames( counts ) <- paste(c(rep(control,n1),rep(exp,n2)),c(1:n1,1:n2),sep="")
#head(counts)
#dim( counts )
#colSums( counts )         ## Library Sizes
#colSums( counts ) / 1e06     ##Library Sizes in millions of reads
#table( rowSums( counts ) )[ 1:30 ]  # Number of genes with low counts
group <- c(rep(control, n1) , rep(exp, n2))
bcv<-0.2
cds <- DGEList( counts , group = group )
#names( cds )
#head(cds$counts)  # original count matrix
#cds$samples  # contains a summary of your samples
#sum( cds$all.zeros ) # How many genes have 0 counts across all samples
#cds # or type the name of the object
# keep only those genes that have at least 1 read per million in at least 2 samples
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 2, ]
#dim( cds )
cds <- calcNormFactors( cds )
#cds$samples
# effective library sizes
#cds$samples$lib.size * cds$samples$norm.factors
# To view the plot immediately
#plotMDS( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )
#crea pdf
#pdf( "/biodata/ArabidopsisVirus/smallRNA/t-test/MDS_CaMV_Mock.pdf" , width = 7 , height = 7 )
#dev.off()
# Estimating Dispersions
#cds <- estimateCommonDisp( cds )
#names( cds )
# The estimate
#cds$common.dispersion

#cds <- estimateTagwiseDisp( cds , prior.df= 10 )
#cds <- estimateTagwiseDisp( cds  )
#names( cds )
#summary( cds$tagwise.dispersion )

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
de.tgw <- exactTest( cds , pair = c( control , exp ),dispersion=bcv^2 )
#de.tgw <- exactTest( cds , pair = c( control , exp ),dispersion=1e-06 )
#de.poi <- exactTest( cds , dispersion = 1e-06 , pair = c( control , exp ) )
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
#head( resultsTbl.tgw)
colnames( resultsTbl.tgw ) <- c( "logFC" , "logConc" , "pVal.Tgw" , "adj.pVal.Tgw" )
write.table( resultsTbl.tgw , file = file_out , sep="\t", row.names = TRUE )




#wh.rows.tgw <- match( rownames( resultsTbl.tgw ) , rownames( cds$counts ) )
# Change column names to be specific to the analysis, logConc and logFC are the same in both.
# Below provides the info to re-order the count matrix to be in line with the order of the results.
#colnames( resultsTbl.poi ) <- c( "logFC" , "logConc" , "pVal.Tgw" , "adj.pVal.Tgw" )
#head( wh.rows.tgw )
# Tagwise Results
#combResults.tgw <- cbind( resultsTbl.tgw ,"Tgw.Disp" = cds$tagwise.dispersion[ wh.rows.tgw ] ,"UpDown.Tgw" = decideTestsDGE( de.tgw , p.value = pval )[ wh.rows.tgw ] ,cds$counts[ wh.rows.tgw , ] )
#head( combResults.tgw )
# Common Results
#write.table( combResults.tgw , file = file_out, sep = "\t" , row.names = TRUE )

#write.table( resultsTbl.tgw , file = file_out, sep = "\t" , row.names = TRUE )


