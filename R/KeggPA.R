#' @name kegg.list
#' @title Kegg Human Annotations Dataframe
#' @descriptiion Obtained from kegg FTP, this dataframe has human annotations for pathway analysis for Entrez IDs and human gene symbols
#' @docType data
#' @usage kegg.list
#' @format First column as Pathway ID, second column as entrez ID in format of "hsa:XXXXX", 3rd column are human gene symbols and 4th column is name of the pathway.
#' @source http://www.genome.jp/kegg/
#' @author Mehran Karimzadeh
NULL



#' Kegg Human Pathway Analysis
#' Based on Fisher's on tail exact test, overrepresentated pathways are found.
#' @author Mehran Karimzadeh
#' @param hgnc a vector of Human Gene symbol IDs or entrez gene IDs in format of "hsa:XXXX"
#' @param project_name a character stating name of the project used in writing the results file in the directory
#' @param kegg.list is a dataframe in format of Kegg Pathway ID, entrez IDs, Gene Symbol IDs and Name of Pathway obtained from kegg directory. By default December 2013 version will be used
#' @return A dataframe of pathway IDs, pathway names, mapped genes, p-values and FDRs.
#' @details Based on assumptions of hypergeometric test, we omit all the genes that are not in kegg database. The hypergeometric p-value is calcualted based on total genes in annotated in kegg database.
#' @export
#' @examples
#' keggPA(c("BRCA1","BRCA2","TP53","XRCC5"),project_name="Example_please_delete")
keggPA <- function(hgnc,project_name, kegg.list=NULL){
  if(is.null(kegg.list)){
    kegg.list <- KeggPA::kegg.list
  }
  if(any(hgnc%in%kegg.list[,2])){
    Z=2
  }else if(any(hgnc%in%kegg.list[,3])){
    Z=3
  }
  KEGG.IDS <- hgnc
  hgnc <- intersect(hgnc,kegg.list[,Z])
  pathways=as.character(levels(as.factor(kegg.list$Name.Pathway)))
  df <- t(sapply(pathways,function(x){
    pathways.genes <- unique(as.character(kegg.list[which(kegg.list[,4]==x),Z]))
    our.genes <- unique(as.character(unique(hgnc)))
    all.kegg <- unique(as.character(kegg.list[,Z]))
    A=length(intersect(our.genes,pathways.genes))###Shared Pathway Query
    C=length(setdiff(our.genes,pathways.genes)) ###Difference query pathway
    B=length(setdiff(pathways.genes,our.genes)) ###
    D=length(all.kegg)-(A+B+C)
    mat <- matrix(c(A,C,B,D),nrow=2)
    colnames(mat) <-  c("Genelist","Genome")
    rownames(mat) <- c("Pathway","Not Pathway")
    p.value <- fisher.test(mat,alternative="greater")[[1]]
    pathway.name <- unique(as.character(kegg.list[which(kegg.list[,4]==x),1]))
    mapped.genes <- paste(unique(intersect(our.genes,pathways.genes)),collapse="|")
    N=A
    total=length(pathways.genes)
    total.queried <- length(unique(KEGG.IDS))
    total.mapped.to.kegg <- length(our.genes)
    result <- c(x,pathway.name,p.value,mapped.genes,N,total,total.queried,total.mapped.to.kegg)
    return(result)
  }))
  df <- as.data.frame(df)
  colnames(df) <- c("Pathway.Name","Pathway.ID","p.value","Mapped.IDs","Number.Mapped.IDs",
                    "Total.Pathway.Genes","Total.genes.queried","Total.genes.Mapped.to.kegg")
  df$BH.FDR <- p.adjust(as.numeric(as.character(df$p.value)),method="BH")	
  df <- df[order(df$p.value,decreasing=F),]
  df <- df[order(df$BH.FDR,decreasing=F),]
  time <- paste(Sys.time())
  time <- unlist(strsplit(time," ",T))[1]
  
  write.table(df,file=paste(project_name,time,"KEGG_Pathway_Analysis.csv"),row.names=F,quote=F,sep="\t")
  return(df)
}		
