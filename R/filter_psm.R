#' filter Peptide
#'
#' @param combined_peptide Combined_peptide file from fragpipe
#' @param combined_modified_peptide combined_modified_peptide file from fragpipe, lbf data no need
#' @param fastafile Fasta files for known proteins
#' @param peptide_quant fragpipe outputs a quantitative file of peptide, for label-free data please use the msstats.csv file
#' @param datatype tmt or lbf
#' @param label Tags are used to filter known proteins, for Swiss-Prot proteins can use 'sp', It is also possible to use species names e.g. 'HUMAN', 'MOUSE', as long as it is possible to distinguish between classical and non-classical proteins.
#'
#' @return Quantification of peptide intensity after filtration
#' @export
#'
#' @examples
#'
#' #TMT data
#'
#'test<-filter_psm('/scratch/lb4489/project/MS_data/test_fdr/tmt_313/combined_peptide.tsv',
#'                 '/scratch/lb4489/project/MS_data/test_fdr/tmt_313/combined_modified_peptide.tsv',
#'                 '/scratch/lb4489/bioindex/uniprot_human.fasta',
#'                 '/scratch/lb4489/project/MS_data/test_fdr/tmt_313/tmt-report/abundance_peptide_MD.tsv',
#'                 'tmt',
#'                 'sp')
#'
#' #label-free data
#'
#'test2<-filter_psm('/scratch/lb4489/project/MS_data/test_fdr/lbf_3_14/combined_peptide.tsv',
#'                  '/scratch/lb4489/project/MS_data/test_fdr/lbf_3_14/combined_modified_peptide.tsv',
#'                  '/scratch/lb4489/bioindex/uniprot_human.fasta',
#'                  '/scratch/lb4489/project/MS_data/test_fdr/lbf_3_14/MSstats.csv',
#'                  'lbf',
#'                  'sp')
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

  nonadd<-nonadd[grep(label,nonadd$Protein,invert =T),]

  list1<-nonadd$Peptide.Sequence[grep(label,nonadd$Mapped.Proteins)]

  nonadd<-nonadd[grep(label,nonadd$Mapped.Proteins,invert =T),]

  if(nrow(nonadd)>1)  {

  peplist<-nonadd$Peptide.Sequence
  peplist<-unique(peplist)

  peplist<-data.frame(seq=c(peplist),
                      match=c(NA))

  pattern2 <- Biostrings::readAAStringSet(fastafile)

  print('Filter all possible single amino acid polymorphism products')

  for (i in 1:nrow(peplist)){

    pattern1 <- Biostrings::AAString(peplist$seq[i])

    match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)

    peplist$match[i]<-sum(match)
  }

  list2<-peplist$seq[peplist$match>0]

  quant<-data.table::fread(peptide_quant,sep = '\t',header = T)

  filter<-unique(c(list1,list2))
  quant<-quant[!(quant$Peptide %in% filter),]

  return(quant) } else{ print('No products found')}

  }

  else{

    pep<-data.table::fread(combined_peptide)
    modpep<-data.table::fread(combined_modified_peptide)

    nonadd<-data.frame(Peptide.Sequence=c(pep$`Peptide Sequence`,modpep$`Peptide Sequence`),
                       Protein=c(pep$Protein,modpep$Protein),
                       Mapped.Proteins=c(pep$`Mapped Proteins`,modpep$`Mapped Proteins`))

    nonadd<-nonadd[grep(label,nonadd$Protein,invert =T),]

    list1<-nonadd$Peptide[grep(label,nonadd$Mapped.Proteins)]

    nonadd<-nonadd[grep(label,nonadd$Mapped.Proteins,invert =T),]

    if(nrow(nonadd)>1) {

    peplist<-nonadd$Peptide
    peplist<-unique(peplist)

    peplist<-data.frame(seq=c(peplist),
                        match=c(NA))

    pattern2 <- Biostrings::readAAStringSet(fastafile)

    print('Filter all possible single amino acid polymorphism products')

    for (i in 1:nrow(peplist)){

      pattern1 <- Biostrings::AAString(peplist$seq[i])

      match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)

      peplist$match[i]<-sum(match)
    }

    list2<-peplist$seq[peplist$match>0]

    quant<-data.table::fread(peptide_quant)

    filter<-unique(c(list1,list2))

    quant<-quant[!(gsub("[^A-Z]", "", quant$PeptideSequence) %in% filter),]

    return(quant) } else{
      print('No products found')
    }
  }
}

#' Filter diann output results
#'
#' @param peptide peptide file from fragpipe
#' @param fastafile Fasta files for known proteins
#' @param report diann output
#' @param label Tags are used to filter known proteins, for Swiss-Prot proteins can use 'sp', It is also possible to use species names e.g. 'HUMAN', 'MOUSE', as long as it is possible to distinguish between classical and non-classical proteins.
#'
#' @return A filtered diann quantitative result, which can be used for subsequent analysis using the diann package.
#' @export
#'
#' @examples
#'
#' test3<-filter_dia_psm('/scratch/lb4489/project/van_cell/cety_regi/fragpipe/peptide.tsv',
#'                      '/scratch/lb4489/bioindex/uniprot_mouse.fasta',
#'                      '/scratch/lb4489/project/van_cell/cety_regi/fragpipe/diann-output/report.tsv',
#'                      'sp')
#'
filter_dia_psm<- function(peptide,fastafile,report,label) {
  nonadd<-data.table::fread(peptide)
  nonadd<-nonadd[-grep(label,nonadd$Protein),]

  list1<-nonadd$Peptide[grep(label,nonadd$`Mapped Proteins`)]

  nonadd<-nonadd[-grep(label,nonadd$`Mapped Proteins`),]

  peplist<-nonadd$Peptide
  peplist<-unique(peplist)

  peplist<-data.frame(seq=c(peplist),
                      match=c(NA))

  pattern2 <- Biostrings::readAAStringSet(fastafile)

  print('Filter all possible single amino acid polymorphism products')

  for (i in 1:nrow(peplist)){

    pattern1 <- Biostrings::AAString(peplist$seq[i])

    match<-Biostrings::vcountPattern(pattern1, pattern2, max.mismatch=1)

    peplist$match[i]<-sum(match)
  }

  list2<-peplist$seq[peplist$match>0]

  quant<-data.table::fread(report)

  filter<-unique(c(list1,list2))

  quant<-quant[!(quant$Stripped.Sequence %in% filter),]

  return(quant)
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


#' Conversion of screened peptide quantification files to MSstats objects for estimating protein intensities. This function only works with label-free data.
#'
#' @param quant Filtered peptide intensity, filter_psm output result
#' @param ...
#'
#' @return A list contains two objects, list[[1]] is an MSstats object that can be used for quality control and subsequent analyses, list[[2]] is the estimated protein intensity obtained using the MSstats algorithm
#' @export
#'
#' @seealso \code{\link{MSstats::dataProcess}}
#'
#' @examples
#'
#' lbf<-build_msstats_lbf(test2)
#'
build_msstats_lbf <- function(quant,...) {

  quant<-as.data.frame(quant)
  quant$BioReplicate[is.na(quant$BioReplicate)]<-1

  print('Converting quantitative results to MSstatsTMT objects can take a long time depending on the number of peptides and sample size')

  processedData <- MSstats::dataProcess(quant, ...) #特别慢这个function

  protein<-MSstats::quantification(processedData)

  list<-list()
  list[[1]]<-processedData
  list[[2]]<-protein

  return(list)
}

#' Use this function to merge data from different fractions
#'
#' @param quant filter_dia_psm output
#' @param meta fp-manifest file
#' @param merge T or F, set T to merge different fractions, if merged a meta file is required.
#'
#' @return Merged protein expression matrix
#' @export
#'
#'
#' @examples
#'
#' dia<-dia_quant(test3,'/scratch/lb4489/project/van_cell/cety_regi/file.fp-manifest',merge=F)
#'
dia_quant <- function (quant,meta,merge) {

  protein.groups <- diann::diann_maxlfq(quant[quant$Q.Value <= 0.01 &quant$PG.Q.Value <= 0.01,],
                                        group.header="Protein.Group", id.header = "Precursor.Id",
                                        quantity.header = "Precursor.Normalised")

  if (!merge) {
    return(protein.groups)
  } else{

  meta<-read.table(meta)
  meta<-meta[meta$V3=='DIA',]

  exp<-as.data.frame(array(NA,c(nrow(protein.groups),length(unique(meta$V2)))))

  colnames(exp)<-unique(meta$V2)
  rownames(exp)<-rownames(protein.groups)

  for (i in 1:ncol(exp)){
    file<-meta$V1[meta$V2==colnames(exp)[i]]

    signs<-protein.groups[,colnames(protein.groups) %in% file]
    signs[is.na(signs)]<-0

    if(length(file)>1) {

      exp[,i]<-rowSums(signs)

    } else{

      exp[,i]<-signs

    }
  }

  return(exp) }

}





