#' Generate gtf files for Ribo-seq analysis
#'
#' @param annotate_files Annotation file, such as gencode.v44.annotation.gff3
#' @param assemble_transcriptome GffCompare or gffcompare results
#'
#' @return Adjusted assembled transcriptome gtf files so that they can be merged with gtf files for Ribo-seq analysis
#' @export
#'
#' @examples
#'
#' test<-generate_annotations(annotate_files='gencode.v44.annotation.gff3',
#'                             assemble_transcriptome='1fpkm.annotated.gtf')
#'
#' write.table(test[,-10],file = 'test.gtf',sep = '\t',quote = F,col.names = F,row.names = F)
#'
#'
generate_annotations <- function(annotate_files, assemble_transcriptome){

  ref <- as.data.frame(rtracklayer::import.gff(annotate_files,format = 'gff'))

  # ref<-ref[ref$type=='gene',]
  #
  # ref<-ref[,c('gene_id','gene_type','gene_name')]

  new <- as.data.frame(rtracklayer::import.gff(assemble_transcriptome,format = 'gff'))

  new<-new[new$strand!='*',]

  new<-new[!(new$transcript_id %in% new$transcript_id[new$class_code=='=']),]

  cache<-as.data.frame(array(NA,c(1,9)))

  print('Start adjusting the annotation file, if there are a large number of assembled transcripts this may take more time')

  library(progress)
  pb <- progress_bar$new(total = length(unique(new$transcript_id)))

  for (i in unique(new$transcript_id)){

    pb$tick()

    inf<-new[new$transcript_id==i,]


    if(!is.na(inf$ref_gene_id[1])) {

      anno<-as.data.frame(array(NA,c(nrow(inf),9)))
      anno$V1<-c(inf$seqnames[1])
      anno$V2<-c('STRINGTIE')
      anno$V3<-c(as.character(inf$type))
      anno$V4<-c(inf$start)
      anno$V5<-c(inf$end)
      anno$V6<-c('.')
      anno$V7<-c(inf$strand[1])
      anno$V8<-c('.')

      transv9=c(paste0('gene_id "',inf$ref_gene_id[1],'";',
                       'transcript_id "',inf$transcript_id[1],'";',
                       'gene_type "',inf$ref_genetype[1],'";',
                       'gene_name "',inf$gene_name[1],'";',
                       'transcript_type "',inf$class_code[1],'";',
                       'transcript_name "',inf$transcript_id[1],'";'))

      exonv9=c(paste0('gene_id "',inf$ref_gene_id[1],'";',
                      'transcript_id "',inf$transcript_id[1],'";',
                      'gene_type "',inf$ref_genetype[1],'";',
                      'gene_name "',inf$gene_name[1],'";',
                      'transcript_type "',inf$class_code[1],'";',
                      'transcript_name "',inf$transcript_id[1],'";',
                      'exon_number "',inf$exon_number,'";',
                      'exon_id "',paste0(inf$transcript_id,'_',inf$exon_numbe),'";'))
      exonv9<-exonv9[-1]

      anno$V9<-c(transv9,exonv9)

    } else{

      anno<-as.data.frame(array(NA,c(nrow(inf)+1,9)))
      anno$V1<-c(inf$seqnames[1])
      anno$V2<-c('STRINGTIE')
      anno$V3<-c('gene',as.character(inf$type))
      anno$V4<-c(inf$start[1],inf$start)
      anno$V5<-c(inf$end[1],inf$end)
      anno$V6<-c('.')
      anno$V7<-c(inf$strand[1])
      anno$V8<-c('.')

      genev9=c(paste0('gene_id "',inf$gene_id[1],'";',
                      'gene_type "',inf$class_code[1],'";',
                      'gene_name "',inf$gene_id[1],'";'))

      transv9=c(paste0('gene_id "',inf$gene_id[1],'";',
                       'transcript_id "',inf$transcript_id[1],'";',
                       'gene_type "',inf$class_code[1],'";',
                       'gene_name "',inf$gene_id[1],'";',
                       'transcript_type "',inf$class_code[1],'";',
                       'transcript_name "',inf$transcript_id[1],'";'))

      exonv9=c(paste0('gene_id "',inf$gene_id[1],'";',
                      'transcript_id "',inf$transcript_id[1],'";',
                      'gene_type "',inf$class_code[1],'";',
                      'gene_name "',inf$gene_id[1],'";',
                      'transcript_type "',inf$class_code[1],'";',
                      'transcript_name "',inf$transcript_id[1],'";',
                      'exon_number "',inf$exon_number,'";',
                      'exon_id "',paste0(inf$transcript_id,'_',inf$exon_numbe),'";'))
      exonv9<-exonv9[-1]

      anno$V9<-c(genev9,transv9,exonv9)

    }

    cache<-rbind(cache,anno)

  }

cache<-cache[-1,]

f <-function(x) unlist(strsplit(x['V9'],';'))[1]
cache$gene <-apply(cache,1,f)

cache2<-cache[1,]


for (i in unique(cache$gene)){
  inf<-cache[cache$gene==i,]

  if(sum(inf$V3 == 'gene')>1) {
    gene<-inf[1,]
    inf<-inf[inf$V3!='gene',]
    gene$V4<-min(inf$V4)
    gene$V5<-max(inf$V5)
    adj<-rbind(gene,inf)
    cache2<-rbind(cache2,adj)
  } else{
    cache2<-rbind(cache2,inf)
  }

}

cache2<-cache2[-1,]

return(cache2)

}

