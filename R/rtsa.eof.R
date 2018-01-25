# title         : Empirical Orthogonal Function analysis of raster time series
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
# 
#' @title EOF (Empirical Orthogonal Functions analysis)
#' 
#' @description This function conducts an Empirical Orthogonal Function analysis (EOF) via a covariance matrix 
#' (cov4gappy function) using "sinkr" package especially designed to handle gappy raster time series
#' 
#' @param rasterts Input raster time series as \code{\linkS4class{RasterStackTS}} or \code{\linkS4class{RasterBrickTS}} object.
#' @param rastermask Either a \code{\linkS4class{RasterLayer}} or "compute". Raster layer to use as a mask. When "compute" 
#' is set raster mask is computed to remove all pixels with incomplete time series.
#' @param nu Numeric. Defines the number of EOFs to return. Defaults to return the full set of EOFs.
#' @param gapfill Character. Defines the algorithm to be used to interpolate pixels with incomplete temporal profiles. Accepts argument supported as method in function \code{\link[rtsa]{rtsa.gapfill}}.
#' @param cores Integer. Defines the number of CPU to be used for multicore processing. Default to "1" core for singlecore processing. Applies only to the masking step.
#' @param centered Logical. If TRUE center the input data before EOF computation. Argument is passed through to function \code{\link[sinkr]{eof}}.
#' @param scaled Logical. If TRUE scale the input data before EOF computation. Argument is passed through to function \code{\link[sinkr]{eof}}.
#' @param ... Additional arguments to be passed through to function \code{\link[sinkr]{eof}}.
#' 
#' @return Object of class \code{\linkS4class{EOFstack}} containing the following components:
#' \tabular{rll}{
#' \tab \code{eof.modes} \tab EOF modes as \code{\linkS4class{RasterBrick}} object\cr
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
#' @references
#' Bjoernsson, H. and Venegas, S.A. (1997). "A manual for EOF and SVD 
#' analyses of climate data", McGill University, CCGCR Report No. 97-1, 
#' Montreal, Quebec, 52pp.
#'
#' Marc, T.H., Losch, M., Wenzel, M., Schroeter, J. (2013). 
#' On the Sensitivity of Field Reconstruction and Prediction Using 
#' Empirical Orthogonal Functions Derived from Gappy Data. Journal of Climate, 
#' 26, 9194-9205. \href{http://dx.doi.org/10.6084/m9.figshare.732650}{pdf}
#' 
#' @keywords EOF PCA SVD time series analysis
#' 
#' @seealso \code{\link[sinkr]{eof}}, \code{\link[rtsa]{rtsa.scaleEOF}}, \code{\link[rtsa]{rtsa.gapfill}}
#' 
#' @examples
#' \dontrun{
#' ## create raster time series using the 'pacificSST' data from 'remote' package
#' require(remote)
#' 
#' data(pacificSST)
#' pacificSST[which(getValues(pacificSST == 0))] <- NA # set NA values
#' # create rts object
#' rasterts <- rts(pacificSST, seq(as.Date('1982-01-15'), as.Date('2010-12-15'), 'months'))
#' 
#' ## generate raster mask
#' raster_mask <- pacificSST[[1]] # create raster mask
#' names(raster_mask) <- "mask"
#' values(raster_mask) <- 1 # set raster mask values
#' raster_mask[which(is.na(getValues(pacificSST[[1]])))] <- 0 # set raster mask values
#'
#' ## compute EOF
#' # compute the first 10 EOFs
#' eof_result <- rtsa.eof(rasterts=rasterts, nu=10)
#' # recursively compute the first 10 EOFs using raster mask
#' eof_result_recursive <- rtsa.eof(rasterts=rasterts, rastermask=raster_mask, nu=10, recursive=T)
#' # compute the first 10 EOFs applying centering, scaling 
#' # and raster mask computing before eof computation
#' eof_res_masked <- rtsa.eof(rasterts=rasterts, rastermask="compute", nu=10, centered=T, scaled=T)
#' }
#' 
#' @import raster
#' @import rts
#' @importFrom sinkr eof
#' @importFrom xts as.xts
#' @importFrom methods new
#' @importFrom stats time
#' @export

# define rsta.eof function
rtsa.eof <- function(rasterts, rastermask=NULL, nu=NULL, gapfill="none", centered=TRUE, scaled=FALSE, method="svds", recursive=FALSE, cores=1L, verbose=FALSE){
  
  # require the 'sinkr' package
  if(!requireNamespace("sinkr", quietly = TRUE))
    stop("Package 'sinkr' is required to run this function.\nYou can install from GitHub repository using the commands:\nlibrary(devtools)\ninstall_github('marchtaylor/sinkr')")

  # check if input file is an object of class 'RasterStackTS', 'RasterBrickTS'
  if(!(class(rasterts) %in% c("RasterStackTS", "RasterBrickTS")))
    stop("'rasterts' argument must be an object of class 'RasterStackTS' or 'RasterBrickTS'.\nUse 'rts()' function to generate 'rasterts' input")
  
  # set number of EOFs to return
  if(is.null(nu)){
    nu <- as.integer(nlayers(rasterts))
    warning(paste("Number of EOFs to return not set using the 'nu' argument.\nIt is set by default to number of input layers: ", as.integer(rasterts@raster@data@nlayers), sep=""))
  } else {
    if(!(class(nu) %in% c("integer", "numeric"))){
      stop("'nu' argument must be numeric")
    } else {
      if(nu > as.integer(nlayers(rasterts)) | nu < 2){
        nu <- as.integer(nlayers(rasterts))
        warning(paste("Number of EOFs to return is higher than the number of input layers.\nIt is set by default to number of input layers: ", as.integer(rasterts@raster@data@nlayers), sep=""))
      }
    }
  }
  if(verbose){
    message("Number of EOFs to compute is: ", nu)
  }
  
  # set raster mask
  # check if 'rastermask' argument is present
  if(!(is.null(rastermask))){
    if(class(rastermask) %in% c("character")){
      if(rastermask %in% c("compute")){
        # computer raster mask from raster time series (pixels)
        if(verbose){
          message("Mask raster will be computed from input raster time series\nPixel temporal profiles with missing data will be masked")
        }
        # create matrix with raster time series values
        matrice_full <- as.matrix(getValues(rasterts))
        
        # define function to mask incomplete pixel temporal profiles
        generateMask <- function(x){
          m <- as.integer(!anyNA(x))
          return(m)
        }
        
        # generate raster mask
        if(cores>1){
          if(!requireNamespace("parallel", quietly = TRUE)){
            warning("Package 'parallel' is required to run this function.\nGoing on using one core")
            generated_mask <- as.vector(apply(X=matrice_full, FUN=generateMask, MARGIN=1))
          } else {
            # check that machine has available cores
            if(cores > as.integer(parallel::detectCores())){
              cores <- as.integer(parallel::detectCores())
            }
            cl <- parallel::makePSOCKcluster(cores)
            parallel::clusterExport(cl, varlist=c("generateMask"), envir=environment())
            generated_mask <- parallel::parApply(cl=cl, X=matrice_full, FUN=generateMask, MARGIN=1)
            parallel::stopCluster(cl)
          }
        } else {
          generated_mask <- as.vector(apply(X=matrice_full, FUN=generateMask, MARGIN=1))
        }
        # set not NaN pixel list
        na_index_mask <- as.vector(which(generated_mask==1))
        # remove imported matrix and clean workspace
        rm(matrice_full)
        gc(verbose=FALSE)
      }
    } else {
      if(!(class(rastermask) %in% c("RasterLayer"))){
        warning("'rastermask' argument is not an object of class 'raster'.\nNo raster mask will be used to clip the input raster time series")
        na_index_mask <- as.vector(1:ncell(rasterts))
      } else {
        # compare mask raster extent with the input raster time series
        mask_correspondence <- as.logical(compareRaster(rasterts@raster, rastermask, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, rotation=TRUE, values=FALSE, stopiffalse=FALSE))
        if(!(mask_correspondence)){
          warning("'rastermask' extent does not correspond to 'rasterts'.\nNo mask will be used to clip the input raster time series")
          na_index_mask <- as.vector(1:ncell(rasterts))
        } else {
          # check mask raster values
          if(length(which(getValues(rastermask)==1))<1){
            warning("'rastermask' does not contain valid pixels for masking purpose (mask pixel values equal to 1).\nNo mask will be used to clip the input raster time series")
            na_index_mask <- as.vector(1:ncell(rasterts))
          } else {
            
            # import raster time series applying raster mask
            if(verbose){
              message("Mask raster is used to mask input raster time series")
            }
            na_index_mask <- as.vector(which(getValues(rastermask)==1))
          }
        }
      }
    }
  } else {
    na_index_mask <- as.vector(1:ncell(rasterts))
    if(verbose){
      message("Mask raster not supplied. Using the full raster time series")
      }
  }
  
  # import raster time series values
  if(verbose){
    message("Importing dataset")
  }
  matrice <- as.matrix(rasterts[na_index_mask])
  if(verbose){
    message("Masked raster time series has ", ncol(matrice), " pixels and ", nrow(matrice), " temporal observations")
  }
  ### matrice object should have time over rows and pixels over columns (contrary of the function help)
  
  # check if there are NAs in the masked dataset
  na_check <- anyNA(matrice)
  if(na_check){
    if(gapfill == "none"){
      warning("Raster time series still contain NA values after masking.\nChange input raster mask or consider the use of the available options\nfor 'gapfill' argument to select a gap-filling method before the EOF computation.\nGoing on however using gappy input raster time series for the EOF computation")
    }
  } else {
    if(verbose){
      message("Raster time series does not contain NA values after masking: OK")
    }
  }

  # perform gap-filling
  if(na_check){
    if(gapfill != "none"){
      # check if gapfill argument is set
      if(!(gapfill %in% c("none", "dineof", "linear", "spline", "stine")))
        stop("'gapfill' argument must be one of the following options: 'none', 'dineof', 'linear', 'spline', 'stine'")
      # generate rastermask layer for gap-filling procedure
      rastermask_gapfill <- raster(rasterts@raster[[1]])
      values(rastermask_gapfill) <- 0
      rastermask_gapfill[na_index_mask] <- 1
      names(rastermask_gapfill) <- "mask"
      # perform gap-filling
      if(verbose){
        message("Gap-filling using '", gapfill, "' method will be applied to masked raster time series")
      }
      rasterts <- rtsa.gapfill(x=rasterts, rastermask=rastermask_gapfill, method=gapfill)
    }
  } else {
    if(verbose){
      warning("Raster time series does not contain NA values after masking: OK.\nGap-filling will be not performed before the EOF computation") 
    }
  }
  
  # check if number of pixel is higher than number of observations
  if(ncol(matrice) <= nrow(matrice)){
    warning("Number of temporal observations is higher than number of masked pixel to be processed.\nThis may generate incorrect results. Going on anyway")
  }
  
  # compute EOF using 'sinkr' package
  if(verbose){
    message(paste(c("EOF analysis started at: "), Sys.time(), sep=""))
  }
  ptm <- proc.time()
  eofresult <- sinkr::eof(matrice, centered=centered, scaled=scaled, method=method, recursive=recursive, nu=nu)
  if(verbose){
    message(paste(c("EOF computation ended at: "), Sys.time(), sep=""))
    message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))
  }
  
  # create output object

  explained_variance <- as.vector(eofresult$Lambda[1:nu]/eofresult$tot.var*100)
  #lambda <- as.vector(eofresult$Lambda[1:nu])
  total_variance <- as.numeric(eofresult$tot.var)
  
  # EOF Expansion Coefficient
  eof.ec <- data.frame(eofresult$A[,1:nu])
  # change column names
  colnamesec <- rep("a", length(eof.ec))
  for(l in 1:length(eof.ec)){
    colnamesec[l] <- paste(c("EC_"), sprintf("%003d", l), sep="")
  }
  names(eof.ec) <- colnamesec
  eof.ec <- as.xts(eof.ec, time(rasterts@time))
  
  # create output raster
  eof_dataset <- brick(rasterts@raster[[1:nu]])
  eof_dataset[na_index_mask] <- eofresult$u ### check if this can be optimized by creating a matrix with the same dimension of the raster data slot
  # generate EOF band names
  eof_names <- as.character(rep(0, nu))
  for(l in 1:nu){
    eof_names[l] <- as.character(paste("EOF mode", sprintf("%003d", l), sep=" "))
  }
  # set band names
  names(eof_dataset) <- validNames(eof_names)
  
  # set final mask
  void_raster <- raster(rasterts@raster[[1]])
  names(void_raster) <- validNames("mask")
  values(void_raster) <- 0
  void_raster[na_index_mask] <- 1
  
  # assemble results in a object of class 'EOFstack'
  eofreturn <- new("EOFstack")
  eofreturn@mask <- void_raster
  eofreturn@eof.modes <- eof_dataset
  eofreturn@expansion_coefficients <- eof.ec
  eofreturn@total_variance <- total_variance
  eofreturn@explained_variance <- explained_variance
  
  # create raster for center values
  if(centered){
    eof.centered <- raster(rasterts@raster[[1]])
    eof.centered[na_index_mask] <- eofresult$F1_center
    eofreturn@center <- eof.centered
  }
  
  # create raster for scaled values
  if(scaled){
    eof.scaled <- raster(rasterts@raster[[1]])
    eof.scaled[na_index_mask] <- eofresult$F1_scale
    eofreturn@scale <- eof.scaled
  }
  
  # return function result
  return(eofreturn)
  
  # stop cluster
  on.exit(stopCluster(cl))
}
