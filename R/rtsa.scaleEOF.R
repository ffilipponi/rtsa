# title         : Scale EOF modes
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
# 
#' @title Scale EOF modes
#' 
#' @description This function rescale resulting spatial modes generated from Empirical Orthogonal Function analysis (EOF)
#' in the range between -1 and 1
#' 
#' @param x Input \code{\link{EOFstack-class}} object generated using \code{\link[rtsa]{rtsa.eof}}
#' @param cut Numeric. Defines the percentile to be used to scale each EOFs. Default is set to 1 (1st and 99th percentiles)
#' @param method Character. Method used to scale EOF modes. Only default \code{"percentile"} is supported
#' 
#' @return Updates the 'eof' slot in the object of class \code{\link{EOFstack-class}} containing the following components:
#' \tabular{lll}{
#' \tab \code{eof} \tab Scaled EOF modes as \code{\linkS4class{RasterBrick}} object\cr
#' \tab \code{expansion_coefficients} \tab EOF Expansion Coefficients (EC) as \code{\linkS4class{xts}} object\cr
#' \tab \code{total_variance} \tab Numeric. Total variance of input raster time series\cr
#' \tab \code{explained_variance} \tab Numeric vector. Percentage of variance explained by each EOF mode with respect to the total variance of input raster time series\cr
#' \tab \code{center} \tab Center values from each pixel temporal profile as \code{\linkS4class{RasterLayer}} object (only computed if \code{centered = TRUE})\cr
#' \tab \code{scale} \tab Scale values from each pixel temporal profile as \code{\linkS4class{RasterLayer}} object (only computed if \code{scaled = TRUE})
#' }
#' 
#' @details 
#' 
#' @author Federico Filipponi
#' 
#' @keywords EOF PCA SVD time series analysis
#' 
#' @seealso \code{\link[rtsa]{rtsa.eof}}, \code{\link[sinkr]{eof}}
#' 
#' @examples
#' \dontrun{
#' ## Scale EOF modes by clipping the distribution between 1st and 99th percentile
#' scaled_eof <- rtsa.scaleEOF(x=eof_result, cut=1)
#' ## Scale EOF modes by clipping the distribution between 5th and 95th percentile
#' scaled_eof <- rtsa.scaleEOF(x=eof_result, cut=5)
#' }
#' 
#' @import raster rts
#' 
#' @export

rtsa.scaleEOF <- function(x, cut=1, method="percentile"){
  
  # check if input data is aobject of class 'EOFstack'
  if(!(class(x) %in% c("EOFstack")))
    stop("'x' argument must be an object of class 'EOFstack' generated using function rtsa.eof()")
  
  if(!(class(x@eof) %in% c("RasterStack", "RasterBrick")))
    stop("'x@eof' argument must be an object of class 'RasterStack' or 'RasterBrick'")
  
  # check if 'cut' argument is in the range 1-99
  cut <- as.integer(round(cut))
  if(cut<1 | cut>99){
    stop("'cut' argument must be in the range 1-99")
  }
  
  # define the number of EOF layers to be processed
  nu <- nlayers(x@eof)
  
  # scale EOF modes
  if(method=="percentile"){
    for(l in 1:nu){ ### valuta se possibile ottimizzare la funzione
      mini <- as.double(quantile(x@eof[[l]], probs=c(0+(cut/100)), na.rm=TRUE))
      massi <- as.double(quantile(x@eof[[l]], probs=c(1-(cut/100)), na.rm=TRUE))
      x@eof[[l]] <- as.vector(ifelse(getValues(x@eof[[l]])>=massi, 1, ifelse(getValues(x@eof[[l]])<= mini, -1, (-1+(((getValues(x@eof[[l]])-mini)/(massi-mini))*2)))))
    }
  }

  # return result
  return(x)
}
