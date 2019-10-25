############################################ NanoPrepareG ############################################



#' @title NanoPrepareG
#' @description Organize MinION/GridION sequencing summary and FASTQ
#' @param DataSummary Path to sequencing summary 
#' @param DataFastq Path to passed FASTQ
#' @details FASTQ are found recursively
#' @return Object of class list
#' @examples
#' #do not run
#' DataSummary<-'path/to/sequencing_summary'
#' DataFastq<-'path/to/fastq_pass'
#' List<-NanoPrepareG(DataSummary,DataFastq)



NanoPrepareG<-function(DataSummary,DataFastq) {

  List<-list()
    
  FastqFiles<-list.files(DataFastq,pattern=".fastq",full.names=TRUE, recursive=TRUE)
  message('Passed FASTQ: ', length(FastqFiles))
    

  Read_Table_Summary<-function(File) {

    Table<-read.table(File,header=FALSE,sep="\t",skip=1)
    RealativeTimeToAdd<-(as.numeric(Table[,11])+as.numeric(Table[,13]))
    Read_Id<-as.character(Table[,3])
    Channel<-as.numeric(Table[,5])
    Mux<-as.numeric(Table[,6])    
    Length<-as.numeric(Table[,14])
    Qscore<-as.numeric(Table[,15])
    Table<-cbind(Read_Id,Channel,Mux,RealativeTimeToAdd,Length,Qscore)
    return(Table)
  }

  message('Reading ', file.path(DataSummary), ' ...')

  SummaryTable<-Read_Table_Summary(DataSummary)

  colnames(SummaryTable)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read","Quality")
  
  List[['fastq']] <-FastqFiles
  List[['summary']]<-SummaryTable

  message("Done")

  return(List)
  
}

