
#' filter Peptide
#'
#' @param combined_peptide
#' @param fastafile
#' @param peptide_quant
#' @param datatype
#'
#' @return Filtered Peptide list for protein intensity estimation
#' @export
#'
#' @examples
#'
#'
#'
filter_psm <- function(combined_peptide,fastafile,peptide_quant,datatype) {
  require('Biostrings')
  require('utils')

  if(datatype=='tmt') {

  nonadd<-utils::read.table(combined_peptide,sep = '\t',header = T)
  nonadd<-nonadd[-grep('HUMAN',nonadd$Protein),]
  nonadd<-nonadd[-grep('sp',nonadd$Protein),]

  list1<-nonadd$Peptide.Sequence[-grep('HUMAN',nonadd$Mapped.Proteins)]

  peplist<-nonadd$Peptide.Sequence
  peplist<-unique(peplist)

  peplist<-data.frame(seq=c(peplist),
                      match=c(NA))

  pattern2 <- Biostrings::readAAStringSet(fastafile)

  for (i in 1:nrow(peplist)){

    pattern1 <- Biostrings::AAString(peplist$seq[i])

    match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)

    peplist$match[i]<-sum(match)
  }

  list2<-peplist$seq[peplist$match>0]

  quant<-utils::read.table(peptide_quant,sep = '\t',header = T)

  filter<-unique(c(list1,list2))
  quant<-quant[!(quant$Peptide %in% filter),]

  return(quant)

  }

  else{

    nonadd<-utils::read.table(combined_peptide,sep = '\t',header = T)
    nonadd<-nonadd[-grep('HUMAN',nonadd$Protein),]
    nonadd<-nonadd[-grep('sp',nonadd$Protein),]

    list1<-nonadd$Peptide.Sequence[-grep('HUMAN',nonadd$Mapped.Proteins)]

    peplist<-nonadd$Peptide.Sequence
    peplist<-unique(peplist)

    peplist<-data.frame(seq=c(peplist),
                        match=c(NA))

    pattern2 <- Biostrings::readAAStringSet(fastafile)

    for (i in 1:nrow(peplist)){

      pattern1 <- Biostrings::AAString(peplist$seq[i])

      match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)

      peplist$match[i]<-sum(match)
    }

    list2<-peplist$seq[peplist$match>0]

    quant<-utils::read.table(peptide_quant,sep = ',',header = T)

    filter<-unique(c(list1,list2))

    quant<-quant[!(gsub("[^A-Z]", "", quant$PeptideSequence) %in% filter),]

    return(quant)
  }
}


#这个代码目前是一个问题，想看看fragpipe团队这么说吧


