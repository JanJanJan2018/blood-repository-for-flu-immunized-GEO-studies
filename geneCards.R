
# These functions grab the genes and the gene summaries from genecards.org
# and some calculate the mean and median fold change values across 
# samples of treatment/control or diseased/healthy etc.

# To extract the genes:
# - find25genes(protein) will grab the 25 genes associated with the protein from web
# - getProteinGenes(protein) will print the genes associated with the protein
# - getSummaries(gene, protein) will grab the gene protein summaries from web
# -getGeneSummaries(protein) will print the gene summary of protein gene
# 

# Many of the flu genes originating file have been downloaded, because when knitr
# was used the rvest web scrape package wouldn't work and error out 

# This file is 'bloodSummaries.csv' to manually extract the gene summaries from 
#12,845 genes.
bloodFluSummaries <- read.csv('bloodSummaries.csv', sep=',',
                              header=TRUE, na.strings=c('',' ','NA'),
                              stringsAsFactors = FALSE)

# The getMeanMedian(gene) will get the mean and median values for all instances of the gene in two samples being compared for pathogenesis/disease/control/treatment/etc

# These are the table names for the functions read in from the specific csv files to use in the functions:
  
FluImz_7day <- read.csv('FluImz_7dayDF.csv',sep=',',header=TRUE,
                        na.strings=c('',' ','NA'),
                        stringsAsFactors=FALSE,row.names=1)
FluImz_1day <- read.csv('FluImz_1dayDF.csv',sep=',',header=TRUE,
                        na.strings=c('',' ','NA'),
                        stringsAsFactors=FALSE,row.names=1)
healthy_noFluImz <- read.csv('healthy_noFluImz.csv',sep=',',header=TRUE,
                             na.strings=c('',' ','NA'),
                             stringsAsFactors=FALSE)
UL <- read.csv('UL_funtionDF.csv',sep=',',header=TRUE,
               na.strings=c('',' ','NA'),
               stringsAsFactors=FALSE)
nonUL <- read.csv('nonUL_functionDF.csv',sep=',',header=TRUE,
                  na.strings=c('',' ','NA'),
                  stringsAsFactors=FALSE)


# This script can be used by calling
# source('geneCards.R') in other script inside same folder. The table are specific
# to the functions to get the fold change values and in this folder.

library(rvest)
library(lubridate)
library(dplyr)

Gene_Path <- './gene scrapes'

if (dir.exists(Gene_Path)){
  unlink(Gene_Path, recursive=TRUE)
  dir.create(Gene_Path)
} else {
  dir.create(Gene_Path)
}

find25genes <- function(protein){
  
  url <- 'https://www.genecards.org/Search/Keyword?queryString=protein'
  
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','%20',protein)
  
  url <- as.character(url)
  url <- gsub('protein',protein, url)
  
  webpage <- read_html(url,encoding = "UTF-8")
  
  protein_html <- html_nodes(webpage,'.symbol-col a')
  protein1 <- html_text(protein_html)
  
  Protein <- as.data.frame(protein1)
  colnames(Protein) <- 'proteinType'
  Protein$proteinType <- as.character(paste(Protein$proteinType))
  Protein$proteinType <- gsub('\n','',Protein$proteinType)
  
  
  date <- as.data.frame(rep(date(),length(Protein$proteinType)))
  colnames(date) <- 'todaysDate'
  
  protein2 <- gsub('%20','-',protein)
  
  proteinName <- as.data.frame(rep(protein2,length(Protein$proteinType)))
  colnames(proteinName) <- 'proteinSearched'
  
  tableProtein <- cbind(Protein,proteinName,date)
  
  setwd(Gene_Path)
  
  
  write.table(tableProtein, 
              paste(protein2,".csv",sep=''), append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  names <- colnames(tableProtein)
  write.csv(names,paste('tableProteinHeader_',protein2,'.csv',sep=''),row.names=FALSE)
  
  setwd('../')
}


#find25genes('estrogen')


getProteinGenes <- function(protein){
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','-',protein)
  table <- read.csv(paste(Gene_Path,'/',protein,'.csv',sep=''),sep=',',
                    header=F,na.strings=c('',' ','NA'), stringsAsFactors = F)
  header <- read.csv(paste(Gene_Path,'/tableProteinHeader_',protein,'.csv',sep=''),
                     sep=',', header=T, na.strings=c('',' ','NA'), stringsAsFactors = F)
  names <- header$x
  colnames(table) <- names
  fileName <- paste('Top25',protein,'s.csv',sep='')
  write.csv(table, fileName, row.names=FALSE)
  return(list(table,fileName))
}


#getProteinGenes('estrogen')



getSummaries <- function(gene,protein){
  url <- 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE&keywords=protein'
  
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ',',',protein)
  gene <- as.character(gene)
  gene <- tolower(gene)
  
  url <- as.character(url)
  url <- gsub('GENE',gene,url)
  url <- gsub('protein',protein, url)
  
  webpage <- read_html(url,encoding = "UTF-8")
  
  Entrez_html <- html_nodes(webpage, '.gc-section-header+ .gc-subsection p')
  Entrez <- html_text(Entrez_html) 
  
  GeneCards_html <- html_nodes(webpage, '.gc-subsection-header+ p')
  GeneCards <- html_text(GeneCards_html) 
  
  UniProt_html <- html_nodes(webpage, '#summaries li:nth-child(1) div')
  UniProtKB <- html_text(UniProt_html) 
  
  Entrez0 <- ifelse(length(Entrez)==0, 'no summary',as.character(paste(Entrez)))
  Entrez1 <- as.data.frame(Entrez0)
  colnames(Entrez1) <- 'EntrezSummary'
  
  GeneCards0 <- ifelse(length(GeneCards)==0,'no summary',
                       as.character(paste(GeneCards)))
  GeneCards1 <- as.data.frame(GeneCards0)
  colnames(GeneCards1) <- 'GeneCardsSummary'
  
  UniProtKB0 <- ifelse(length(UniProtKB)==0,'no summary',
                       as.character(paste(UniProtKB)))
  UniProtKB1 <- as.data.frame(UniProtKB0)
  colnames(UniProtKB1) <- 'UniProtKB_Summary'
  
  Entrez1$EntrezSummary <- as.character(paste(Entrez1$EntrezSummary))
  Entrez1$EntrezSummary <- gsub('\n','',Entrez1$EntrezSummary)
  
  GeneCards1$GeneCardsSummary <- as.character(paste(GeneCards1$GeneCardsSummary))
  GeneCards1$GeneCardsSummary <- gsub('\n','',GeneCards1$GeneCardsSummary)
  
  UniProtKB1$UniProtKB_Summary <- as.character(paste(UniProtKB1$UniProtKB_Summary))
  UniProtKB1$UniProtKB_Summary <- gsub('\n','',UniProtKB1$UniProtKB_Summary)
  
  date <- as.data.frame(rep(date(),length(Entrez1$EntrezSummary)))
  colnames(date) <- 'todaysDate'
  
  protein2 <- gsub(',','-',protein)
  
  proteinName <- as.data.frame(rep(protein2,length(Entrez1$EntrezSummary)))
  colnames(proteinName) <- 'proteinSearched'
  
  gene <- as.data.frame(rep(toupper(gene),length(Entrez1$EntrezSummary)))
  colnames(gene) <- 'gene'
  
  tableProtein <- cbind(proteinName,gene,Entrez1,GeneCards1,UniProtKB1,date)
  
  setwd(Gene_Path)
  
  
  write.table(tableProtein, 
              paste(protein2,"summary.csv",sep=''), append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  names <- colnames(tableProtein)
  write.csv(names,paste('geneHeader_summary_',protein2,'.csv',sep=''),row.names=FALSE)
  
  setwd('../')
  return(gene)
}


#getSummaries('TP53','estrogen')


getGeneSummaries <- function(protein){
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','-',protein)
  
  table <- read.csv(paste(Gene_Path,'/',protein,'summary.csv',sep=''),
                    sep=',',header=F,na.strings=c('',' ','NA'), stringsAsFactors = F)
  
  header <- read.csv(paste(Gene_Path,'/geneHeader_summary_',protein,'.csv',sep=''),
                     sep=',', header=T, na.strings=c('',' ','NA'), stringsAsFactors = F)
  names <- header$x
  colnames(table) <- names
  
  fileName <- paste('proteinGeneSummaries_',protein,'.csv',sep='')
  write.csv(table, fileName, row.names=FALSE)
  return(list(table,fileName))
}


#getGeneSummaries('estrogen')


# These two files are the files that the getMeanMedian('gene') will calculate
# the foldchange mean and median values of all instances of that gene.

UL <- read.csv('UL_funtionDF.csv',sep=',',header=TRUE,
               na.strings=c('',' ','NA'),
               stringsAsFactors=FALSE)
nonUL <- read.csv('nonUL_functionDF.csv',sep=',',header=TRUE,
                  na.strings=c('',' ','NA'),
                  stringsAsFactors=FALSE)




if (dir.exists('./UL and nonUL foldchange tables')){
  unlink('./UL and nonUL foldchange tables', recursive=TRUE)
  dir.create('./UL and nonUL foldchange tables')
} else {
  dir.create('./UL and nonUL foldchange tables')
}


getMeanMedian <- function(gene){
  gene <- as.character(paste(gene))
  gene0_ul <- UL[grep(gene,UL$Gene.Symbol),]
  gene0_nonul <- nonUL[grep(gene,UL$Gene.Symbol),]
  
  sub_ul <- subset(gene0_ul, gene0_ul$Gene.Symbol==gene |
                     gene0_ul$first==gene |
                     gene0_ul$third==gene |
                     gene0_ul$second==gene)
  
  sub_nonul <- subset(gene0_nonul, gene0_nonul$Gene.Symbol==gene|
                        gene0_nonul$first==gene |
                        gene0_nonul$third==gene |
                        gene0_nonul$second==gene)
  
  gene1_UL <- sub_ul[,2:6]
  gene1_nonUL <- sub_nonul[,2:6]
  
  gene1_UL$mean <- apply(gene1_UL,1,mean)
  gene1_UL$median <- apply(gene1_UL,1,median)
  gene1_nonUL$mean <- apply(gene1_nonUL,1,mean)
  gene1_nonUL$median <- apply(gene1_nonUL,1,median)
  
  gene1_UL$FoldChange_mean <- gene1_UL$mean/gene1_nonUL$mean
  gene1_UL$FoldChange_median <- gene1_UL$median/gene1_nonUL$median
  
  geneMeans <- gene1_UL$FoldChange_mean
  
  geneMedians <- gene1_UL$FoldChange_median
  
  
  print('The foldchage of UL means to nonUL means is:')
  print(geneMeans)
  
  print('The foldchage of UL medians to nonUL medians is:')
  print(geneMedians)
  
  colnames(gene1_UL) <- paste(colnames(gene1_UL), '_UL')
  colnames(gene1_nonUL) <- paste(colnames(gene1_nonUL), '_nonUL')
  
  setwd('./UL and nonUL foldchange tables')
  
  write.table(gene1_UL[2:length(gene1_UL$median),], "allGenesUL.csv", append=TRUE, 
              col.names=FALSE, sep=",",
              row.names=TRUE)
  UL_names <- colnames(gene1_UL)
  write.csv(UL_names,'header_UL_names.csv',row.names=FALSE)
  
  write.table(gene1_nonUL[2:length(gene1_nonUL$median),], "allGenesNonUL.csv", append=TRUE, 
              col.names=FALSE, sep=",",
              row.names=TRUE)
  nonUL_names <- colnames(gene1_nonUL)
  
  write.csv(nonUL_names,'header_nonUL_names.csv', row.names=FALSE)
  
  setwd('../')
  
  return(list(gene1_UL,gene1_nonUL))
}


#getMeanMedian("TF")

###################################################################################
# 
# #Lets look at the top 25 genes for 'tumor' and then get the summaries.
# find25genes('tumor')
# 
# #Print out the top 25 genes for the protein searched.
# getProteinGenes('tumor')
# 
# #Grab the gene summaries for a particular gene and the protein interested in.
# getSummaries('TP53','tumor')
# 
# #Print the results for one of the genes and the protein searched.
# getGeneSummaries('tumor')
# 

# Now lets look at how this gene does in UL tissue compared to non-UL tissue as far as fold change goes, with the ration of gene expression in UL/non-UL samples as a median and mean.

# getMeanMedian('TP53')


#Lets use a different protein as double word search in genecards.org, 'hair loss' to be exact.
# find25genes('hair loss')
# 
# getProteinGenes('hair loss')
# 
# getSummaries('GJB2','hair loss')
# 
# getGeneSummaries('hair loss')

#####################################################################################

FluImz_7day <- read.csv('FluImz_7dayDF.csv',sep=',',header=TRUE,
                        na.strings=c('',' ','NA'),
                        stringsAsFactors=FALSE,row.names=1)
FluImz_1day <- read.csv('FluImz_1dayDF.csv',sep=',',header=TRUE,
                        na.strings=c('',' ','NA'),
                        stringsAsFactors=FALSE,row.names=1)
healthy_noFluImz <- read.csv('healthy_noFluImz.csv',sep=',',header=TRUE,
                             na.strings=c('',' ','NA'),
                             stringsAsFactors=FALSE)
UL <- read.csv('UL_funtionDF.csv',sep=',',header=TRUE,
               na.strings=c('',' ','NA'),
               stringsAsFactors=FALSE)
nonUL <- read.csv('nonUL_functionDF.csv',sep=',',header=TRUE,
                  na.strings=c('',' ','NA'),
                  stringsAsFactors=FALSE)


getMeanMedianFlu1 <- function(gene){
  gene <- as.character(paste(gene))
  gene0_ul <- FluImz_1day[grep(gene,FluImz_1day$Gene.Symbol),]
  gene0_nonul <- healthy_noFluImz[grep(gene,healthy_noFluImz$Gene.Symbol),]
  
  sub_ul <- subset(gene0_ul, gene0_ul$Gene.Symbol==gene )
  
  sub_nonul <- subset(gene0_nonul, gene0_nonul$Gene.Symbol==gene)
  
  gene1_UL <- sub_ul[,2:4]
  gene1_nonUL <- sub_nonul[,2:7]
  
  gene1_UL$mean <- apply(gene1_UL,1,mean)
  gene1_UL$median <- apply(gene1_UL,1,median)
  gene1_nonUL$mean <- apply(gene1_nonUL,1,mean)
  gene1_nonUL$median <- apply(gene1_nonUL,1,median)
  
  gene1_UL$FoldChange_mean <- gene1_UL$mean/gene1_nonUL$mean
  gene1_UL$FoldChange_median <- gene1_UL$median/gene1_nonUL$median
  
  geneMeans <- gene1_UL$FoldChange_mean
  
  geneMedians <- gene1_UL$FoldChange_median
  
  
  print('The foldchage of flu immunized after 1 day means to heathy non-immunized flu means is:')
  print(geneMeans)
  
  print('The foldchage of flu immunized after 1 day medians to non-immunized flu medians is:')
  print(geneMedians)
  
  colnames(gene1_UL) <- paste(colnames(gene1_UL), '_flu_1day')
  colnames(gene1_nonUL) <- paste(colnames(gene1_nonUL), '_healthyNonImmz')
  
  
  write.table(gene1_UL[2:length(gene1_UL$median),], "allFluImz1day.csv", append=TRUE, 
              col.names=FALSE, sep=",",
              row.names=TRUE)
  UL_names <- colnames(gene1_UL)
  write.csv(UL_names,'header_allFluImz1day.csv',row.names=FALSE)
  
  write.table(gene1_nonUL[2:length(gene1_nonUL$median),], "allhealthyNonImzFlu.csv", append=TRUE, 
              col.names=FALSE, sep=",",
              row.names=TRUE)
  nonUL_names <- colnames(gene1_nonUL)
  
  write.csv(nonUL_names,'header_allhealthyNonImzFlu.csv', row.names=FALSE)
  
  
  return(list(gene1_UL,gene1_nonUL))
}




getMeanMedianFlu7 <- function(gene){
  gene <- as.character(paste(gene))
  gene0_ul <- FluImz_7day[grep(gene,FluImz_7day$Gene.Symbol),]
  gene0_nonul <- healthy_noFluImz[grep(gene,healthy_noFluImz$Gene.Symbol),]
  
  sub_ul <- subset(gene0_ul, gene0_ul$Gene.Symbol==gene )
  
  sub_nonul <- subset(gene0_nonul, gene0_nonul$Gene.Symbol==gene)
  
  gene1_UL <- sub_ul[,2:4]
  gene1_nonUL <- sub_nonul[,2:7]
  
  gene1_UL$mean <- apply(gene1_UL,1,mean)
  gene1_UL$median <- apply(gene1_UL,1,median)
  gene1_nonUL$mean <- apply(gene1_nonUL,1,mean)
  gene1_nonUL$median <- apply(gene1_nonUL,1,median)
  
  gene1_UL$FoldChange_mean <- gene1_UL$mean/gene1_nonUL$mean
  gene1_UL$FoldChange_median <- gene1_UL$median/gene1_nonUL$median
  
  geneMeans <- gene1_UL$FoldChange_mean
  
  geneMedians <- gene1_UL$FoldChange_median
  
  
  print('The foldchage of flu immunized after 7 days means to heathy non-immunized flu means is:')
  print(geneMeans)
  
  print('The foldchage of flu immunized after 7 days medians to non-immunized flu medians is:')
  print(geneMedians)
  
  colnames(gene1_UL) <- paste(colnames(gene1_UL), '_flu_7day')
  colnames(gene1_nonUL) <- paste(colnames(gene1_nonUL), '_healthyNonImmz')
  
  
  write.table(gene1_UL[2:length(gene1_UL$median),], "allFluImz7day.csv", append=TRUE, 
              col.names=FALSE, sep=",",
              row.names=TRUE)
  UL_names <- colnames(gene1_UL)
  write.csv(UL_names,'header_allFluImz7day.csv',row.names=FALSE)
  
  write.table(gene1_nonUL[2:length(gene1_nonUL$median),], "allhealthyNonImzFlu.csv", append=TRUE, 
              col.names=FALSE, sep=",",
              row.names=TRUE)
  nonUL_names <- colnames(gene1_nonUL)
  
  write.csv(nonUL_names,'header_allhealthyNonImzFlu.csv', row.names=FALSE)
  
  
  return(list(gene1_UL,gene1_nonUL))
}


#####################################################################################


