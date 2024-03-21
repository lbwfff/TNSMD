#' Co-expression analysis using the GSEA algorithm
#'
#' @param array Protein intensity array
#' @param name Name of peptide analysed
#' @param species Only 'human' or 'mouse' is supported.
#'
#' @return A list consists of [1] a clusterProfiler object, [2] a GSEA diagram
#' @export
#'
#' @examples
coexp_gsea<-function(array,name,species) {

  library(clusterProfiler)
  library(ggplot2)

  if (species=='human'){organism='hsa'
  library(org.Hs.eg.db)
  Org = org.Hs.eg.db
  } else {organism='mmu'
  library(org.Mm.eg.db)
  Org = org.Mm.eg.db
  }

  pro<-array

  cor<-pro[grep(name,rownames(pro)),]

  exp<-as.data.frame(pro)
  cor<-as.data.frame(cor)

  exp<-exp[grep('rev',rownames(exp),invert=T),]

  corresult<-cor(t(exp),t(cor),use = "p")
  corresult<-as.data.frame(corresult)
  corresult$what<-c('whatever')

  colnames(corresult)[1]<-c('sel')
  corresult<-corresult[order(corresult$sel,decreasing = T),]
  corresult<-corresult[-1,]

  idmatch <-clusterProfiler::bitr(rownames(corresult), fromType = "UNIPROT",
                 toType = c( "ENTREZID"),
                 OrgDb = Org)
  idmatch<-idmatch[match(rownames(corresult),idmatch$UNIPROT),]
  corresult$ID<-idmatch$ENTREZID

  corresult<-corresult[!is.na(corresult$ID),]

  enrich<-as.numeric(corresult$sel)
  names(enrich) <- as.character(corresult$ID)


  edo2 <- clusterProfiler::gseGO(enrich,OrgDb= Org,ont='ALL',pvalueCutoff = 0.5,eps=0)

  p<-list()

  p[[1]]<-
    GseaVis::gseaNb(object = edo2,
           geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='BP'][1:5],
           subPlot = 2,
           termWidth = 35,
           # legend.position = c(0.8,0.8),
           addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
           pvalX = 0.05,pvalY = 0.05)

  p[[2]]<-
    GseaVis::gseaNb(object = edo2,
           geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='MF'][1:5],
           subPlot = 2,
           termWidth = 35,
           # legend.position = c(0.8,0.8),
           addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
           pvalX = 0.05,pvalY = 0.05)

  p[[3]]<-
    GseaVis::gseaNb(object = edo2,
           geneSetID = edo2@result$ID[edo2@result$ONTOLOGY=='CC'][1:5],
           subPlot = 2,
           termWidth = 35,
           # legend.position = c(0.8,0.8),
           addPval = T,curveCol=c("#76BA99", "#EB4747", "#996699","#08519C", "#A50F15"),
           pvalX = 0.05,pvalY = 0.05)
  list<-list(object=edo2,plot=p)
  return(list)
}

