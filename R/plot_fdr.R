
#' Visualisation of subgroup FDR control
#'
#' @param target target psm file, e.g. 10NIH_M_451_032523_percolator_target_psms.tsv
#' @param decoy decoy psm file, e.g. 10NIH_M_451_032523_percolator_decoy_psms.tsv
#' @param label Label for the classical proteome, 'sp' when using the uniprot proteome
#'
#' @return A list consists of three ggplot images.
#' @export
#'
#' @examples
plot_fdr <- function(target,decoy,label) {

  target<-read.table(target,header = T)
  decoy<-read.table(decoy,header = T)

  plot<-data.frame(
    score=c(target$score,decoy$score),
    group=c(rep('target',nrow(target)),rep('decoy',nrow(decoy)))
  )

  library(ggplot2)

  p<-list()

  p[[1]]<-
    ggplot(plot, aes(score, fill = group, col = I("black"))) +
    geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
    labs(x = 'score', y = "",
         title = "Global") +
    scale_fill_manual(values = c("decoy" = "#FF9900","target" = "#009900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)


  target_clas<-target[grep(label,target$proteinIds),]
  decoy_clas<-decoy[grep(label,decoy$proteinIds),]


  plot<-data.frame(
    score=c(target_clas$score,decoy_clas$score),
    group=c(rep('target',nrow(target_clas)),rep('decoy',nrow(decoy_clas)))
  )

  p[[2]]<-
    ggplot(plot, aes(score, fill = group, col = I("black"))) +
    geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
    labs(x = 'score', y = "",
         title = "Canonical protein") +
    scale_fill_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)

  target_clas<-target[-grep(label,target$proteinIds),]
  decoy_clas<-decoy[-grep(label,decoy$proteinIds),]


  plot<-data.frame(
    score=c(target_clas$score,decoy_clas$score),
    group=c(rep('target',nrow(target_clas)),rep('decoy',nrow(decoy_clas)))
  )

  p[[3]]<-
    ggplot(plot, aes(score, fill = group, col = I("black"))) +
    geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
    labs(x = 'score', y = "",
         title = "sORF encoded peptide") +
    scale_fill_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)

return(p)
}

#' Visualising the distribution of scores for specific peptides
#'
#' @param target target psm file, e.g. 10NIH_M_451_032523_percolator_target_psms.tsv
#' @param decoy decoy psm file, e.g. 10NIH_M_451_032523_percolator_decoy_psms.tsv
#' @param label Label for the classical proteome, 'sp' when using the uniprot proteome
#' @param peptidelist A sequence of peptides, e.g. peptidelist<-c('LIEFCISSMRN','MTMLADHAARQ')
#'
#' @return A list consists of two ggplot images.
#' @export
#'
#' @examples
plot_peptide_fdr<-function(target,decoy,label,peptidelist){

  target<-read.table(target,header = T)
  decoy<-read.table(decoy,header = T)

  target_clas<-target[-grep(label,target$proteinIds),]
  decoy_clas<-decoy[-grep(label,decoy$proteinIds),]

  plot<-data.frame(
    score=c(target_clas$score,decoy_clas$score),
    pep=runif(nrow(target_clas)+nrow(decoy_clas)),
    group=c(rep('target',nrow(target_clas)),rep('decoy',nrow(decoy_clas))),
    seq=c(target_clas$peptide,decoy_clas$peptide)
  )
  plot$score<-as.numeric(plot$score)
  plot$pep<-as.numeric(plot$pep)

  plot$seq <- gsub("[^A-Z]", "", plot$seq)

  library(ggrepel)

  p<-list()

  p[[1]]<-
    ggplot(plot, aes(score, fill = group, col = I("black"))) +
    geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
    labs(x = 'score', y = "",
         title = "sORF encoded peptide") +
    scale_fill_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)



  p[[2]]<-
    ggplot(data = plot,aes(x=score,y=pep, colour=group))+
    geom_point(shape=21,size=4)+
    scale_colour_manual(values = c("target" = "#009900", "decoy" = "#FF9900")) +
    geom_text_repel(data = plot[plot$seq %in% peptidelist,],aes(label = seq),color = "black",
                    size = 3,segment.color = "black", show.legend = FALSE,hjust = -1.5 )+
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1.2)),
          axis.text = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 0))+
    theme(aspect.ratio=1)+
    labs(x = 'score', y = "")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  return(p)
}


