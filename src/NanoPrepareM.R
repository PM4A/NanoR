############################################ NanoPrepareM ############################################


#' @title NanoPrepareM
#' @description Organize MinION/GridION basecalled FAST5
#' @param DataPass Path to passed FAST5
#' @param DataFail Path to failed FAST5 [NA]
#' @param DataSkip Path to skipped FAST5 [NA]
#' @param MultiRead Logical. If TRUE, assume multi-read FAST5 [FALSE]
#' @details FAST5 are found recursively
#' @return Object of class list
#' @examples
#' #do not run
#' DataPass<-"/path/to/fast5_pass"
#' DataFail<-"/path/to/fast5_fail" #can be omitted
#' #single-read .fast5 files
#' List<-NanoPrepareM(DataPass,DataFail,Label=Label)
#' #multi-read .fast5 files
#' List<-NanoPrepareM(DataPass,DataFail, Label=Label,MultiRead=TRUE)



NanoPrepareM<-function(DataPass,DataFail=NA,DataSkip=NA, MultiRead=FALSE) { #store informations about paths and numbers of failed/skipped FAST5
  
  List<-list()
  
  PassFiles<-list.files(DataPass, full.names=TRUE, recursive = TRUE, pattern=".fast5")
  message('Passed FAST5: ', length(PassFiles))
  
  if (is.na(DataFail)) {
    
    FailFilesLength<-0
    message('Failed FAST5: ', FailFilesLength)
    
  }
  
  else {
    
    if (MultiRead==FALSE) {
      
      FailFilesLength<-length(list.files(DataFail, recursive = TRUE, pattern=".fast5"))
      message('Failed FAST5: ', FailFilesLength)
      
    }
    
    else {
      
      library(rhdf5)
      
      FailFiles<-list.files(DataFail, recursive = TRUE, pattern=".fast5", full.names=TRUE)
      message('Failed FAST5: ', length(FailFiles))
      
      FailFilesOrdered<-FailFiles[order(as.numeric(gsub("[^0-9]+", "", FailFiles)))]
      
      First<-FailFilesOrdered[1]
      FileOpen<-H5Fopen(First)
      howmany<-nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)
      
      Last<-FailFilesOrdered[length(FailFilesOrdered)]
      FileOpen<-H5Fopen(Last)
      howmany_<-nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)
      
      FailFilesLength<-(howmany*(length(FailFiles)-1) + howmany_)
    }
    
  }
  
  if (is.na(DataSkip)) {
    
    SkipFilesLength<-0
    message('Skipped FAST5: ', SkipFilesLength)
    
  }
  
  else {
    
    if (MultiRead==FALSE) {
      
      SkipFilesLength<-length(list.files(DataSkip, recursive = TRUE, pattern=".fast5"))
      message('Skipped FAST5: ', SkipFilesLength)
      
    }
    
    else {
      
      library(rhdf5)
      
      SkipFiles<-list.files(DataSkip, recursive = TRUE, pattern=".fast5", full.names=TRUE)
      message('Skipped FAST5: ', length(SkipFiles))
      
      SkipFilesOrdered<-SkipFiles[order(as.numeric(gsub("[^0-9]+", "", SkipFiles)))]
      
      First<-SkipFilesOrdered[1]
      FileOpen<-H5Fopen(First)
      howmany<-nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)
      
      Last<-SkipFilesOrdered[length(SkipFilesOrdered)]
      FileOpen<-H5Fopen(Last)
      howmany_<-nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)
      
      SkipFilesLength<-(howmany*(length(SkipFiles)-1) + howmany_)
      
    }
    
  }
  
  List[['fast5']]<-PassFiles
  List[['failed']]<-FailFilesLength
  List[['skipped']]<-SkipFilesLength
  List[['multiread']]<-MultiRead
  
  message("Done")  
  
  return(List)
  
}