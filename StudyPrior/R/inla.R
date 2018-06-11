# Check if inla is available

check.inla <- function(){
  if(!requireNamespace("INLA", quietly = TRUE)){
    stop("This functionality relies on INLA. Please check it is available in your R library.\n See http://www.r-inla.org/download or run install.inla() \n");
    }
}


#' Install R-INLA from http://www.r-inla.org/download
#'
#' @param testing install the testing version. Default (FALSE) installs the stable version.
#'
#' @return TRUE if INLA can be successfully loaded after installation.

install.inla <- function(testing=FALSE){
  INLA.url = if(testing) "https://inla.r-inla-download.org/R/testing" else "https://inla.r-inla-download.org/R/stable"
  
  print(paste0("Installing INLA with: install.packages('INLA', repos=c(getOption('repos'), INLA=",INLA.url,",dep=TRUE)"))
  utils::install.packages("INLA", repos=c(getOption("repos"), INLA=INLA.url), dep=TRUE)
  return(requireNamespace("INLA"))
}


