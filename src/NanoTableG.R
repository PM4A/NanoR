
############################################ NanoTableG ############################################


#' @title NanoTableG
#' @description Generate a metadata table
#' @param NanoGList Object of class list from NanoPrepareG
#' @param DataOut Output folder
#' @param GCC  Logical. If TRUE, compute GC content from FASTQ [FALSE]
#' @return Object of class matrix
#' @examples
#' #do not run
#' #assume List is the output from NanoPrepareG
#' DataOut <- "/path/to/output"
#' Table<-NanoTableG(List,DataOut)
#' #include GC content
#' Table<-NanoTableM(List,DataOut,GCC=TRUE) 


NanoTableG<-function(NanoGList,DataOut,GCC=FALSE) {

  Directory<-file.path(DataOut)
  dir.create(Directory,showWarnings=FALSE, recursive=TRUE)

  TableInDirectory<-list.files(Directory,pattern="metadata.txt")

  if(length(TableInDirectory) != 0) {

    stop("Cannot use a directory that already contains other results")
  
  }

  Read_Id<-as.character(NanoGList$summary[,1])
  Channel<-as.numeric(NanoGList$summary[,2])
  Mux<-NanoGList$summary[,3]
  Length<-as.numeric(NanoGList$summary[,5])
  Qscore<-as.numeric(NanoGList$summary[,6])
  Relative_Time<-as.numeric(NanoGList$summary[,4])

  if (GCC == TRUE) {
      
    library(ShortRead)

    message("Calculating GC content...") ## very hard to speed up Fastq parsing in R ... :( 
    # 1 hour for 550 fastq file, 8000 sequences each. Definitely slow.
    
    
    FastqFilesPath<-NanoGList[[1]]
    
    GCC<-function(seq) {      
      GC<-sum(gregexpr('[GgCc]',seq)[[1]] > 0)/nchar(seq) #faster than using library
      return(GC)
    }    
    
    Gc_Con<-function(Element) {
      fqFile<-FastqFile(Element)
      Fastq <- tryCatch({
        readFastq(fqFile)},
        error = function(cond) {
          return(NULL)},
        warning = function(cond) {
          message(cond)
          return(NULL)}
      )
      if (is.null(Fastq)) {
        warning("Ill-formatted FASTQ : ", Element, " .Skipped")
      }
      else {
        CharRead<-as.character(sread(Fastq))
        close(fqFile)
        GC<-lapply(CharRead,GCC)
      }
      return(GC)
    }

    List<-lapply(FastqFilesPath,Gc_Con)      
    GC_Content<-unlist(List)
    GCL<-length(GC_Content)
    NT<-nrow(NanoGList[[2]])
    
    if (GCL != NT) { ##lack of some FASTQ . Get GC content from what we have.
      LackOfReads<-NT-GCL
      GC_To_Add<-rep(NA, LackOfReads)
      GC_Content<-c(GC_Content,GC_To_Add)
    }
    
    Table_Tot<-cbind(Read_Id,Channel,Mux,Relative_Time,Length,Qscore,GC_Content)
    colnames(Table_Tot)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read","Quality", "GC content")
    write.table(Table_Tot, file.path(Directory, 'metadata.txt'), sep="\t", quote=FALSE, col.names=T, row.names=FALSE) 
  
  }
  
  else {
        
    GC_Content<-rep("GC_Content",nrow(NanoGList[[2]]))
    Table_Tot<-cbind(Read_Id,Channel,Mux,Relative_Time,Length,Qscore,GC_Content)
    colnames(Table_Tot)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read", "Quality", "GC content")
    write.table(Table_Tot, file.path(Directory, 'metadata.txt'), sep="\t", quote=FALSE, col.names=T, row.names=FALSE) 
  }

  message("Done")
  return(Table_Tot)
  
}

