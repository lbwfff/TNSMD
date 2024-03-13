#' Construction of a proteogenomics index based on ribo-seq results
#'
#' @param file Ribo-seq analysis output file, When source is set to other, you need to provide three columns including 'genename', 'orftype' and 'pepseq'. csv file
#' @param source ribocode, ribotish or other, Can be set to other to index result files or public database files that are not ribocode or ribotish outputs
#' @param length Peptide length, used to filter long proteins in the index, set to 0 when there is no need to filter by length
#' @param label fasta head tag, e.g. "sPep", not 'sp' or 'tr' is recommended.
#' @param OrganismName OrganismName is the scientific name of the organism for the UniProtKB entry, defaulting to "Homo sapiens".
#' @param OrganismIdentifier OrganismIdentifier is the unique identifier of the source organism, assigned by NCBI, defaulting to "9696".
#'
#' @return A Biostrings object can be written as a fasta file
#' @export
#'
#' @examples generate_index(file='ribocode.txt',source='ribocode',length=100,label='sPep')
#'
#'
generate_index<- function(file,source,length,label,
                          OrganismName,OrganismIdentifier) {

  add_numbers <- function(vec) {
    result <- character(length(vec))
    counts <- list()

    for (i in seq_along(vec)) {
      element <- vec[i]
      if (element %in% names(counts)) {
        counts[[element]] <- counts[[element]] + 1
      } else {
        counts[[element]] <- 1
      }
      result[i] <- paste(element, counts[[element]], sep = " ")
    }

    return(result)
  }

  if(source=='ribocode'){

  ribo<-read.table(file,sep = '\t',header = T)
  ribo<-ribo[!duplicated(ribo$ORF_ID),]
  ribo<-ribo[ribo$ORF_type!='annotated',]

  if(length==0) {
    ribo<-ribo
  } else{
    ribo<-ribo[ribo$ORF_length/3 <= length,]
  }

  #构建肽索引了，重点是fasta的head
  orftype<-ifelse(ribo$ORF_type=='dORF','dORF',ifelse(
    ribo$ORF_type=='internal','internal ORF',ifelse(
      ribo$ORF_type=='novel','novel ORF',ifelse(
        ribo$ORF_type=='Overlap_dORF','Overlap_dORF',ifelse(
          ribo$ORF_type=='Overlap_uORF','Overlap_uORF', 'uORF')))))

  ribo$name1<-paste0(label,'|',label,'_',1:nrow(ribo),'|',label,'_',1:nrow(ribo))

  ribo$name2<-paste0(ribo$gene_name,' ',orftype)

  ribo$name2<-add_numbers(ribo$name2)
  if (missing(OrganismName)) {
  ribo$name3<-paste0('OS=','Homo sapiens',' OX=','9606',' GN=',ribo$gene_name,' PE=6')
  } else{
  ribo$name3<-paste0('OS=',OrganismName,' OX=',OrganismIdentifier,' GN=',ribo$gene_name,' PE=6')
  }

 index<-data.frame(name=c(paste0(ribo$name1," ",ribo$name2," ",ribo$name3)),
                   sequence=c(ribo$AAseq))

  }

  if(source=='ribotish'){
    ribo<-data.table::fread(file,sep = '\t',header = T)
    ribo<-ribo[ribo$TisType!='Annotated',]

    if(length==0) {
      ribo<-ribo
    } else{
      ribo<-ribo[ribo$AALen <= length,]
    }

    orftype<-ifelse(ribo$TisType %in% c("3'UTR"),'dORF',ifelse(
      ribo$TisType %in% c("5'UTR"),'uORF',ifelse(
        ribo$TisType %in% c("3'UTR:CDSFrameOverlap"),'Overlap_dORF',ifelse(
          ribo$TisType %in% c("5'UTR:CDSFrameOverlap"),'Overlap_uORF',ifelse(
            ribo$TisType %in% c("Novel"),'novel ORF', 'other ORF')))))

    ribo$name1<-paste0(label,'|',label,'_',1:nrow(ribo),'|',label,'_',1:nrow(ribo))

    ribo$name2<-paste0(ribo$Symbol,' ',orftype)

    ribo$name2<-add_numbers(ribo$name2)
    if (missing(OrganismName)) {
      ribo$name3<-paste0('OS=','Homo sapiens',' OX=','9606',' GN=',ribo$Symbol,' PE=6')
    } else{
      ribo$name3<-paste0('OS=',OrganismName,' OX=',OrganismIdentifier,' GN=',ribo$Symbol,' PE=6')
    }

    index<-data.frame(name=c(paste0(ribo$name1," ",ribo$name2," ",ribo$name3)),
                      sequence=c(ribo$Seq))
  }

  if(source=='other'){
  ribo<-read.csv(file)
  ribo<-ribo[,c('genename','orftype','pepseq')]
  ribo$length<-nchar(ribo$pepseq)

  if(length==0) {
    ribo<-ribo
  } else{
    ribo<-ribo[ribo$length <= length,]
  }

  ribo$name1<-paste0(label,'|',label,'_',1:nrow(ribo),'|',label,'_',1:nrow(ribo))

  ribo$name2<-paste0(ribo$genename,' ',ribo$orftype)

  ribo$name2<-add_numbers(ribo$name2)
  if (missing(OrganismName)) {
    ribo$name3<-paste0('OS=','Homo sapiens',' OX=','9606',' GN=',ribo$genename,' PE=6')
  } else{
    ribo$name3<-paste0('OS=',OrganismName,' OX=',OrganismIdentifier,' GN=',ribo$genename,' PE=6')
  }

  index<-data.frame(name=c(paste0(ribo$name1," ",ribo$name2," ",ribo$name3)),
                    sequence=c(ribo$pepseq))

  }
  return(index)
}

