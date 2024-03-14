#' filter Peptide
#'
#' @param combined_peptide Combined_peptide file from fragpipe
#' @param combined_modified_peptide combined_modified_peptide file from fragpipe
#' @param fastafile Fasta files for known proteins
#' @param peptide_quant fragpipe outputs a quantitative file of peptide
#' @param datatype tmt or lbf
#' @param label Tags are used to filter known proteins, for uniprot proteins use 'sp'
#'
#' @return Quantification of peptide intensity after filtration
#' @export
#'
#' @examples
#'
#'
#'
filter_psm <- function(combined_peptide,combined_modified_peptide,fastafile,
                       peptide_quant,datatype,label) {
  require('Biostrings')
  require('utils')

  if(datatype=='tmt') {

  pep<-data.table::fread(combined_peptide)
  modpep<-data.table::fread(combined_modified_peptide)

  nonadd<-data.frame(Peptide.Sequence=c(pep$`Peptide Sequence`,modpep$`Peptide Sequence`),
                     Protein=c(pep$Protein,modpep$Protein),
                     Mapped.Proteins=c(pep$`Mapped Proteins`,modpep$`Mapped Proteins`))

  nonadd<-nonadd[!duplicated(nonadd$Peptide.Sequence),]

  nonadd<-nonadd[-grep(label,nonadd$Protein),]

  list1<-nonadd$Peptide.Sequence[grep(label,nonadd$Mapped.Proteins)]

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

  quant<-data.table::fread(peptide_quant,sep = '\t',header = T)

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




#' Conversion of screened peptide quantification files to MSstatsTMT objects for estimating protein intensities. This function only works with tmt data.
#'
#' @param quant Filtered peptide intensity, filter_psm output result
#' @param annotation experiment_annotation.tsv file from fragpipe
#'
#' @return  A list contains two objects, list[[1]] is an MSstatsTMT object that can be used for quality control and subsequent analyses, list[[2]] is the estimated protein intensity obtained using the MSstatsTMT algorithm
#' @export
#'
#' @examples This function only works with tmt data.
#'
#'
build_msstats <- function (quant,annotation) {

  data<-as.data.frame(quant)
  data<-data[,c(3,4,10:ncol(data))]

  data<-tidyr::gather(data, key = "sample",
               value = "Intensity", -c('ProteinID', 'Peptide'),
               na.rm = FALSE, convert = FALSE, factor_key = FALSE)

  data<-as.data.frame(data)
  data<-data[!is.na(data$Intensity),]

  annation<-read.table(annotation,header = T)
  annation$run<-rep(paste0('run',annation$plex))
  annation$condition[annation$condition=='Bridge']<-c('Norm')
  annation$replicate<-c(1)
  annation$replicate[annation$condition=='Norm']<-c('Norm')

  input<-data.frame(ProteinName=c(data$ProteinID),PeptideSequence=c(data$Peptide),
                    Charge=c(NA),PSM=c(NA),Mixture=c('Mixture1'),TechRepMixture=c(1),
                    Run=c(NA),Channel=c(NA),BioReplicate=c(NA),
                    Condition=c(NA),Intensity=c(data$Intensity),sample=c(data$sample))

  match<-annation[match(input$sample,annation$sample),]

  input$Run<-match$run
  input$Channel<-match$channel
  input$BioReplicate<-match$replicate
  input$Condition<-match$condition
  input$Intensity<-2^input$Intensity

  input$Charge<-c(1)
  input$PSM<-input$PeptideSequence

  print('Converting quantitative results to MSstatsTMT objects can take a long time depending on the number of peptides and sample size')

  quant.msstats <- MSstatsTMT::proteinSummarization(input,
                                        method="msstats",
                                        global_norm=TRUE,
                                        reference_norm=TRUE,
                                        remove_norm_channel = TRUE,
                                        remove_empty_channel = TRUE)

  protein<-quant.msstats[["ProteinLevelData"]]
  protein$sample<-NA

  for (i in unique(paste0(protein$Run,'_',protein$Channel))){
    sample<-annation$sample[paste0(annation$run,'_',annation$channel)==i]
    protein$sample[paste0(protein$Run,'_',protein$Channel)==i]<-sample
  }

  exparray <- tidyr::spread(protein[,c(5,6,9)], key = sample, value = Abundance)

  list<-list()
  list[[1]]<-quant.msstats
  list[[2]]<-exparray

  return(list)

}


