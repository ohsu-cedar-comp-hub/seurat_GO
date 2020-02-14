args <- commandArgs()

help <- function(){
    cat("runGOforDESeq2.R :
- For the deseq2 output in the pipeline, run GO analysis on significant genes.
- Currently this is compatible with genome assemblies hg19 (Ens75), hg38 (Ens89) and hg38 (Ens90)
- Input DESEq2 table must be in .txt format
- Color options can be hex followed by saturation ex. #FE117A60 or rcolors
- The plot will have the same name as the degFile but with a .pdf extension.\n")
    cat("Usage: \n")
    cat("--degFile  : deseq2 table with log2FoldChange and pvalue [ required ]\n")
    cat("--adjp     : FDR adjusted p-value cutoff                 [ default = 0.01 ]\n")
    cat("--assembly : genome assembly                             [ requires ]\n")
    cat("--FC       : fold change cutoff (not log2 transformed)   [ default = 2 ]\n")
    cat("\n")
    q()
}


degFile = snakemake@input[['degFile']]


assembly <- snakemake@params[['assembly']]

FC <- snakemake@params[['FC']]
adjp <- snakemake@params[['adjp']]


## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args))){
    help()
}

if (identical(adjp,character(0))){
   adjp<-0.01
}else{
    adjp <- as.numeric(adjp)
}

if (identical(FC,character(0))){
    FC <- 1
} else{
    FC <- as.numeric(FC)
}

library(GO.db)
library(topGO)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(GenomicFeatures)
##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)
deg <- read.delim(file=degFile,header=TRUE,row.names=1,sep="\t")
deg <- deg[is.finite(deg$avg_logFC),]
deg <- deg[!grepl('mt-',row.names(deg)),]

##---------load correct Biomart------------------------#
if (assembly == "hg19") {
    organismStr <- "hsapiens"
    geneID2GO <- get(load("/home/groups/CEDAR/anno/biomaRt/hg19.Ens_75.biomaRt.GO.external.geneID2GO.RData"))
    xx <- get(load("/home/groups/CEDAR/anno/biomaRt/GO.db.Term.list.rda"))
}
if (assembly == "mm10.Ens_96") {
  organismStr <- "mm10"
  ### to get to mm10 mappings ensembl 96!
  geneID2GO <- get(load("/home/groups/CEDAR/anno/biomaRt/mm10.Ens_96.biomaRt.external.geneID2GO.RData"))
  xx <- get(load("/home/groups/CEDAR/anno/biomaRt/GO.db.Term.list.rda"))
}

##-----------------------------------Functions--------------------------------------#
runGO <- function(geneList,xx=xx,otype,setName,fname){
  setLength       <- sum(as.numeric(levels(geneList))[geneList]) 
  GOData          <- new("topGOdata", ontology=otype, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)##run topGO
  resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
  x               <- GenTable(GOData, classicFisher=resultFisher, topNodes=length(names(resultFisher@score)))## make go table for all terms
  x               <- data.frame(x)
  pVal            <- data.frame(pval=signif(resultFisher@score, 6)) ## get unrounded pvalue
  x$enrich        <- x$Significant/x$Expected ## calculate enrichment based on what you expect by chance
  x$p.unround     <- pVal[x$GO.ID,"pval"]## put unrounded pvalue in the table
  x$p.adj         <- signif(p.adjust(x$p.unround, method="BH"), 6)## calculate the adjusted pvalue with Benjamini & Hochberg correction
  x$log.p.adj     <- -log10(x$p.adj) ## convert adjusted p value to -log10 for plot magnitude
  x <- x[order(x$GO.ID),]
  write.table(x, file=fname, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE) ## save the table
  ## you can print the tree if you want, but since I keep the list of all of them skip
  
  printGraph(GOData,## make the tree for the go data
             resultFisher,
             firstSigNodes = 5,
             fn.prefix = paste(setName, setLength, otype, sep="_"),
             useInfo = "all",
             pdfSW = TRUE
  )    
  
  return(x)  
}

## function to make barplot of -log10 adjusted pvalues colored by enrichment
drawBarplot <- function(go, ontology, setName, setSize){
  go <- go[!go$p.adj > 0.01,]
  if(nrow(go)>1){
    go$Term <- make.unique(paste(sapply(strsplit(as.character(substring(go$Term,1,50)), "\\,"), `[`, 1)))
    print(setName)
    print(setSize)
    go <- go[with(go, order(p.adj, -enrich)),]
    ## Currently there is a discrepency between xx and x, so we only use Term right now, not Term.full
    go$Term <-factor(paste(go$Term), levels=rev(paste(go$Term))) ## sort table by adjusted p-value
    ptitle <- paste(ontology, setName, sep='_') ## plot title
    ptitle <- gsub("^.*/","",ptitle)
    pfname <- paste(setName,ontology, setSize,"pdf",sep=".")## name of png file
    if(nrow(go) < 20 ){
      toprange <- 1:nrow(go)
    }else{
      toprange <- 1:20
    }
    top <- go[toprange,]
    col <- colorRampPalette(c("white","navy"))(16)   
    pdf(file=paste(Dir, pfname, sep="/"), height=5,width=7)
    print({
      p <- ggplot(top, aes(y=log.p.adj, x=Term, fill=enrich)) + ## ggplot barplot function
        geom_bar(stat="identity",colour="black") +
        ggtitle(ptitle) +
        xlab("") + ylab("-log10(fdr)") +
        scale_fill_gradient(low=col[2], high=col[15], name="enrichment", limits=c(0,ceiling(max(top$enrich))))+
        coord_flip()+
        theme(panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        theme(text = element_text(size=8),
              axis.text.x = element_text(vjust=1,color="black",size=8),
              axis.text.y = element_text(color="black",size=8),
              plot.title=element_text(size=10))   
    })
    dev.off()
  }
}

writeGOTable <- function(geneList,xx=xx,otype,setName, geneSel){
    fname <- paste(Dir, paste(setName, otype, "GO_consolidated.txt", sep="_"), sep="/")
    GOData = new("topGOdata", ontology=otype, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5, geneSel=geneSel)
    resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")
    goresults = GenTable(GOData,classicFisher=resultFisher,topNodes=20)
    goresults$genes <- sapply(goresults$GO.ID, function(x)
        {
          genes<-genesInTerm(GOData, x)
          genes[[1]][genes[[1]] %in% geneSel]
        })
     goresults$genes <-vapply(goresults$genes, paste, collapse=", ", character(1L))
    write.table(goresults, file=fname, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
}



down_file <- snakemake@output[['down']]
up_file <- snakemake@output[['up']]

print("get up genes and make geneList")
up <- deg$p_val_adj < adjp & deg$avg_logFC >= log2(FC)
up <- unique(rownames(deg[up,]))
all <-unique(names(geneID2GO))
up.geneList <-  factor(as.integer(all %in% up))
names(up.geneList) <- all
up.setsize <- sum(as.numeric(levels(up.geneList))[up.geneList])

print("setsize for significant genes") 
up.setsize
adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.txt$|\\.rda$","",degFile)

Dir <- sub("$", "/GOterms", dirname(comparison))
if(!(file.exists(Dir))) {
  dir.create(Dir,FALSE,TRUE)
}

print("make GO table for the up genes")
#################################

if(up.setsize <2){
write.table(NA, file=up_file, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)}else{
go.UP.BP <- runGO(geneList=up.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"upFC",FC, "adjp", adjp, sep="."),fname=up_file)
writeGOTable(geneList=up.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"upFC",FC, "adjp", adjp, sep="."),geneSel=up)
}
print("make the png for the up genes")

if(up.setsize >2){
go.UP.BP <- go.UP.BP[is.finite(go.UP.BP$enrich),]
setName <- paste(substr(comparison,26,nchar(comparison)-4),"barplot","upFC",FC, "adjp", adjp, sep=".")
drawBarplot(go=go.UP.BP,ontology="BP",setName=setName, setSize=up.setsize)
print("get down genes and make geneList")
}

dn <- deg$p_val_adj < adjp & deg$avg_logFC <= -log2(FC)
dn <- unique(rownames(deg[dn,]))
all <-unique(names(geneID2GO))
dn.geneList <-  factor(as.integer(all %in% dn))
names(dn.geneList) <- all
dn.setsize <- sum(as.numeric(levels(dn.geneList))[dn.geneList])
print("setsize for significant genes") 

dn.setsize
print("make GO table for down genes")
if(dn.setsize <2){
write.table(NA, file=down_file, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)}else{
go.DN.BP <- runGO(geneList=dn.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"downFC",FC, "adjp", adjp, sep="."),fname=down_file)
writeGOTable(geneList=dn.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"downFC",FC, "adjp", adjp, sep="."),geneSel=dn)
}

print("make barplot for down genes")
if(dn.setsize >2){
go.DN.BP <- go.DN.BP[is.finite(go.DN.BP$enrich),]
setName <- paste(substr(comparison,26,nchar(comparison)-4),"barplot","downFC",FC, "adjp", adjp, sep=".")
drawBarplot(go=go.DN.BP,ontology="BP",setName=setName, setSize=dn.setsize )}
