# Asthma-biomarkers
To identify biomarkers related to asthma severity 

RMA normalization of microarray data

# Install Packages: 

BiocManager::install("affy") 
BiocManager::install("hgu133a.db") 
BiocManager::install("hgu133acdf")



# Calling Libraries:


library(affy)
library(hgu133a.db)
library(hgu133acdf)


# Set working Directory:

setwd("~/Desktop/asthma_research_file")          # the name of folder which contains ".CELL" raw file

# Data Normalization:

f <-list.files(pattern = ".CEL")
f   # it will list the name of all .CEL files in the folder



j = 0                             # using a loop
i = 0

for (i=1, i<=j,j++ in 1:10)        # n=10 or the n is equal to number of .CELL files to be normalized

{
 raw.data <- ReadAffy(verbose = FALSE, filenames = f[I], cdfname = "hgu133acdf")
 data.rma.norm = rma(raw.data)

}


rma = exprs(data.rma.norm)
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t")




# Annotation:

probes=row.names(rma)
 
ls("package:hgu133a.db")

Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))

Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))

rma=cbind(probes,Symbols,Entrez_IDs,rma)

write.table(rma, file = "annotation.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


