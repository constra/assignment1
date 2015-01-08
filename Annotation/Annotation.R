##Annotation lessons
library("hgu95av2.db")
ids <- c("39730_at", "1635_at", "1674_at", "40504_at", "40202_at")
columns(hgu95av2.db)
keytypes(hgu95av2.db)

library("yeast2.db")
columns(yeast2.db)
keytypes(yeast2.db)
setdiff(keytypes(hgu95av2.db),keytypes(yeast2.db))

head(keys(yeast2.db,keytype="ENZYME"))
head(keys(hgu95av2.db,keytype="ENTREZID"))

##mapping
select(hgu95av2.db,keys=ids,columns="ENTREZID",keytype="PROBEID")

#example
head(keys(yeast2.db,keytype="ENSEMBL"))
gene_names=c("PUF3","PUF4","TOR2")
entrez_id=select(yeast2.db,gene_names,columns="ENTREZID",keytype="GENENAME")
select(yeast2.db,keys=gene_names,columns="ENZYME",keytype="GENENAME")
puf34=select(yeast2.db,keys=gene_names,columns="GO",keytype="GENENAME")

#pass GO term to get GO-id
library(GO.db)
keytypes(GO.db)
select(GO.db,keys=puf34$GO,column="TERM",keytype="GOID")

### Using OrgDB
library(org.Hs.eg.db)
keys=head(keys(org.Hs.eg.db,keytype="ENTREZID"),n=2)
columns=c("PFAM","GO","SYMBOL")
select(org.Hs.eg.db,keys,columns,keytype="ENTREZID")

library(org.Sc.sgd.db)
keytypes(org.Sc.sgd.db)
columns=c("PFAM","ENZYME","GENENAME")
select(org.Sc.sgd.db,gene_names,columns,keytype="GENENAME")
### Exercises
#Exercise 1: Look at the help page for the different columns and keytypes 
#values with: help("SYMBOL"). Now use this information and what we just 
#described to look up the entrez gene, probe id and chromosome for the gene
#symbol "MSX2".
columns=c("ENTREZID","PROBEID","CHR")
select(hgu95av2.db,"MSX2",columns,keytype="SYMBOL")

#Exercise 2: Examine the gene symbols for both the hgu95av2.db and the 
#org.Hs.eg.db packages. Which one has more gene symbols? Which one has more 
#gene symbols that can be mapped to an entrez gene ID? Which object seems to
#contain more information?
symbol_hgu=keys(hgu95av2.db,keytype="SYMBOL")
symbol_org=keys(org.Hs.eg.db,keytype="SYMBOL")
entrez_map_hgu=select(hgu95av2.db,symbol_hgu,"ENTREZID",keytype="SYMBOL")
entrez_map_org=select(org.Hs.eg.db,symbol_org,"ENTREZID",keytype="SYMBOL")
dim(entrez_map_hgu)
length(columns(org.Hs.eg.db)) < length(columns(hgu95av2.db))
#Exercise 3: In the previous exercise we had to use gene symbols as keys. 
#But in the past this kind of behavior has sometimes been inadvisable because 
#some gene symbols are used as the official symbol for more than one gene. 
#To learn if this is still happening take advantage of the fact that entrez gene
#ids are uniquely assigned, and extract all of the gene symbols and their associated
#entrez gene ids from the org.Hs.eg.db package. Then check the symbols for redundancy.
entrezid_hgu=select(hgu95av2.db,symbol_hgu,"ENTREZID","SYMBOL")
length(unique(entrezid_hgu$ENTREZID))
length(unique(entrezid_hgu$SYMBOL))
badids=entrezid_hgu$SYMBOL[duplicated(entrezid_hgu$SYMBOL)]
select(hgu95av2.db,badids,"ENTREZID","SYMBOL")

###using TxDB
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hs=TxDb.Hsapiens.UCSC.hg19.knownGene
keys=head(keys(txdb_hs,keytype="GENEID",n=2))
keytypes(txdb_hs)
columns=c("TXNAME","TXSTART","TXSTRAND")
select(txdb_hs,keys,columns,keytype="GENEID")
transcripts(txdb_hs)
exons(txdb_hs)
cds(txdb_hs)
transcripts(txdb_hs, columns = c("tx_id", "tx_name", "gene_id"))
transcriptsBy(txdb_hs,by="gene")

#Exercise 4: Use the accessors for the TxDb.Hsapiens.UCSC.hg19.knownGene package
#to retrieve the gene id, transcript name and transcript chromosome for all the 
#transcripts. Do this using both the select() method and also using the transcripts() 
#method. What is the difference in the output?
keys=keys(txdb_hs,keytype="GENEID",n=2)
select(txdb_hs,keys,columns=c("TXNAME","TXCHROM"),keytype="GENEID")
transcripts(txdb_hs,columns=c("gene_id","tx_name"))
#Exercise 5: Load the TxDb.Athaliana.BioMart.plantsmart16 package. This package is not from 
#UCSC and it is based on plantsmart. Now use select or one of the range based accessors to 
#look at the gene ids from this TranscriptDb object. How do they compare to what you saw in
#the TxDb.Hsapiens.UCSC.hg19.knownGene package?
library(TxDb.Athaliana.BioMart.plantsmart16)
#
#
#
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
txdb_sc=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
select(txdb_sc,"YPR116W",columns=c("TXNAME","TXCHROM"),keytype="GENEID")
promoters(txdb_sc,upstream=2000,downstream=200)

###using OrganismDB
library(Homo.sapiens)
#library(Saccharomyces.cerevisiae) <- Make one
gd <- list(join1 = c(GO.db = "GOID", org.Sc.sgd.db= "GO"), 
           join2 = c(org.Sc.sgd.db = "ENTREZID",
                     TxDb.Scerevisiae.UCSC.sacCer3.sgdGene = "GENEID"))
makeOrganismPackage(pkgname = "Saccharomyces.cerevisiae", graphData = gd, organism = "Saccharomyces.cerevisiae",
                    version = "1.0.0", maintainer = "Mai Sun<mai.sun@embl.de>", author = "Mai Sun", destDir = ".", license = "Artistic-2.0")


#Exercise 6: Use the Homo.sapiens object to look up the gene symbol, transcript start and chromosome using select(). 
#Then do the same thing using transcripts. You might expect that this call to transcripts will look the same as it 
#did for the TranscriptDb object, but (temporarily) it will not.
keys=keys(txdb_hs,keytype="GENEID",n=2)
select(Homo.sapiens,keys,columns=c("SYMBOL","TXSTART","TXCHROM"),keytype="GENEID")
transcripts(Homo.sapiens,columns=c("SYMBOL","TXSTART","TXCHROM"))
#Exercise 7: Look at the results from call the columns method on the Homo.sapiens object and 
#compare that to what happens when you call columns on the org.Hs.eg.db object and then 
#look at a call to columns on the TxDb.Hsapiens.UCSC.hg19.knownGene object. What is the difference
#between TXSTART and CHRLOC? Which one do you think you should use for transcripts or other genomic information?
columns(Homo.sapiens)
columns(org.Hs.eg.db)
columns(TxDb.Hsapiens.UCSC.hg19.knownGene)
###using AnnotationHub
library(AnnotationHub)
ah=AnnotationHub()
res=ah$goldenpath.hg19.encodeDCC.wgEncodeUwTfbs.wgEncodeUwTfbsMcf7CtcfStdPkRep1.narrowPeak_0.0.1.RData

#Exercise 8: Set the AnnotationHub filter to NULL to clear it out, and then set ip up so that 
#it is extracting data that originated with the UCSC data provider and that is also limited to 
#Homo sapiens and the hg19 genome.
filters(ah)=NULL
columns(ah)
filters(ah)=list(DataProvider="hgdownload.cse.ucsc.edu",Species="Homo sapiens",Genome="hg19")
length(ah)
#Exercise 9: Now that you have basically narrowed things down to the hg19 annotations from UCSC genome
#browser, lets get one of these annotations. Now tab complete your way to the oreganno track and save it into a local variable.
res_oregano=ah$goldenpath.hg19.database.oreganno_0.0.1.RData
###Using biomaRt
library(biomaRt)
head(listMarts)





RNAinteractome=read.csv2(file="DATA/mRNAinteractome",sep="\t",header=T)
binders=as.vector(RNAinteractome[,"name"])
querys=c("AVO1","AVO2","BIT61","LST8","SLM1","SLM2","TOR2","TSC11")
intersect(binders,querys)
library("yeast2.db")
name_ensembl=select(yeast2.db,keys=querys,columns=c("ENSEMBL"),keytype="GENENAME")

library("biomaRt")
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useMart("ensembl",dataset="scerevisiae_gene_ensembl")
listFilters(ensembl)
listAttributes(ensembl)
listAttributes(ensembl)[grep("length",listAttributes(ensembl)[,"description"],ignore.case=T),]
listFilters(ensembl)[grep("ensembl",listFilters(ensembl)[,"description"],ignore.case=T),]

getBM(attributes = "cds_length", filters = "ensembl_gene_id",values = name_ensembl[,"ENSEMBL"], mart = ensembl)

