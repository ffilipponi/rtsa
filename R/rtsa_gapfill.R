# title         : Gap-filling raster time series
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
#
#' @title Raster time series gap-filling
#' 
#' @description This function perform gap-filling of gappy raster time series
#' 
#' @param x Input raster time series as \code{\linkS4class{RasterStackTS}} or \code{\linkS4class{RasterBrickTS}} object
#' @param rastermask A \code{\linkS4class{RasterLayer}} to use as a mask. If not set
#' a raster mask is computed to remove all pixels with less than two values in temporal profiles
#' @param method Character. Defines the algorithm to be used to interpolate pixels with incomplete temporal profiles. 
#' Accepts the following input:
#' \tabular{lll}{
#' \tab \code{"linear"} \tab for linear interpolation in \code{\link[imputeTS]{na.interpolation}} using \code{\link[stats]{approxfun}}\cr
#' \tab \code{"spline"} \tab for spline interpolation in \code{\link[imputeTS]{na.interpolation}} using \code{\link[stats]{splinefun}}\cr
#' \tab \code{"stine"} \tab for stine interpolation in \code{\link[imputeTS]{na.interpolation}} using \code{\link[stinepack]{stinterp}}\cr
#' \tab \code{"dineof"} \tab for dineof interpolation using \code{\link[sinkr]{dineof}}\cr
#' }
#' @param cores Integer. Defines the number of CPU to be used for multicore processing. Default to "1" core for 
#' singlecore processing.
#' @param ... Additional arguments
#' 
#' @return Object of class \code{\linkS4class{RasterBrickTS}} with gap-filled pixels
#' 
#' @details 
#' 
#' @author Federico Filipponi
#' 
#' @references
#' 
#' @keywords time series analysis gap-filling
#' 
#' @seealso \code{\link[imputeTS]{na.interpolation}}, \code{\link[sinkr]{dineof}}, \code{\link[stats]{approxfun}}, \code{\link[stats]{splinefun}}, \code{\link[stinepack]{stinterp}}
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
#' values(raster_mask) <- 1 # set raster mask values
#' raster_mask[which(is.na(getValues(pacificSST[[1]])))] <- 0 # set raster mask values
#' 
#' ## randomly remove values from cells in rts object
#' frac_gaps <- 0.5 # the fraction of data with NaNs
#' temporal_cells <- as.integer(ncell(rasterts) * nlayers(rasterts)) # number of total cells in rts
#' # define random position of cells to be set to NaN
#' na_cells <- sort(unique(sample.int(temporal_cells, (temporal_cells * frac_gaps))))
#' gappy_values <- as.vector(getValues(rasterts)) # extract raster values
#' gappy_values[na_cells] <- NA # set NA to random positions
#' rasterts_gappy <- setValues(rasterts, values=gappy_values) # set NA to pixels
#' 
#' ## perform gap-filling on the gappy dataset
#' 
#' # using linear interpolation
#' rasterts_linear <- rtsa.gapfill(rasterts_gappy, method="linear")
#' 
#' # using spline interpolation and multiple cores
#' rasterts_spline <- rtsa.gapfill(rasterts_gappy, method="spline", cores=4)
#' 
#' # using stine interpolation and raster mask
#' rasterts_stine <- rtsa.gapfill(rasterts_gappy, rastermask=raster_mask, method="stine")
#' 
#' # using dineof interpolation and raster mask
#' rasterts_dineof <- rtsa.gapfill(rasterts_gappy, rastermask=raster_mask, method="dineof")
#' }
#' 
#' @import raster rts sinkr imputeTS
#' @importFrom zoo as.yearmon index
#' @importFrom xts periodicity
#' @importFrom sp coordinates
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' 
#' @export

rtsa.gapfill <- function(x, rastermask=NULL, method, cores=1L, verbose=FALSE){
  
  if(!(class(x) %in% c("RasterStackTS", "RasterBrickTS")))
    stop("'x' argument must be an object of class 'RasterStackTS', 'RasterBrickTS'")
  
  if(!(method %in% c("dineof", "linear", "spline", "stine")))
    stop("'method' argument must be one of the following options: 'dineof', 'linear', 'spline', 'stine'")
  
  #################
  # Process raster mask
  if(!(is.null(rastermask))){
    if(!(class(rastermask) %in% c("RasterLayer"))){
      warning("'rastermask' argument is not an object of class 'raster'.\nNo raster mask will be used to clip the input raster time series")
      na_index_mask <- as.vector(1:as.integer(ncell(x)))
    } else {
      # compare mask raster extent with the input raster time series
      mask_correspondence <- as.logical(compareRaster(x@raster, rastermask, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, rotation=TRUE, values=FALSE, stopiffalse=FALSE))
      if(!(mask_correspondence)){
        warning("'rastermask' extent does not correspond to 'x'.\nNo mask will be used to clip the input raster time series")
        na_index_mask <- as.vector(1:as.integer(ncell(x)))
      } else {
        # check mask raster values
        if(length(which(getValues(rastermask)==1))<1){
          warning("'rastermask' does not contain valid pixels for masking purpose (mask pixel values equal to 1).\nNo mask will be used to clip the input raster time series")
          na_index_mask <- as.vector(1:as.integer(ncell(x)))
        } else {
          
          # import raster time series applying raster mask
          if(verbose){
            message("Mask raster is used to mask input raster time series")
          }
          na_index_mask <- as.vector(which(getValues(rastermask)==1))
        }
      }
    }
  } else {
    if(verbose){
      message("Mask raster not supplied. Masking pixel with less than two temporal observation") ### probabilmente questo check va effettuato anche se la maschera viene fornita ma non maschera tutti i pixel con profili nulli
    }

    # compute mask using pixel that have less than two values in temporal profiles
    # define function to mask incomplete pixel profiles
    generateMask <- function(x){
      m <- as.integer(length(which(is.na(x)))<=(length(x)-2))
      return(m)
    }
    
    # generate raster mask
    if(cores>1){
      if(!requireNamespace("parallel", quietly = TRUE)){
        warning("Package 'parallel' is required to run this function.\nGoing on using one core")
        generated_mask <- as.vector(apply(X=as.matrix(getValues(x)), FUN=generateMask, MARGIN=1))
      } else {
        # check that machine has available cores
        if(cores > as.integer(parallel::detectCores())){
          cores <- as.integer(parallel::detectCores())
        }
        cl <- parallel::makePSOCKcluster(cores)
        parallel::clusterExport(cl, varlist=c("generateMask"), envir=environment())
        generated_mask <- parallel::parApply(cl=cl, as.matrix(getValues(x)), FUN=generateMask, MARGIN=1)
        parallel::stopCluster(cl)
      }
    } else {
      generated_mask <- as.vector(apply(X=as.matrix(getValues(x)), FUN=generateMask, MARGIN=1))
    }
    # set not NaN pixel list
    na_index_mask <- as.vector(which(generated_mask==1))
  }
  
  #################
  # Perform gap-filling using 'dineof' method from package 'sinkr'
  if(method == "dineof"){
    if(!requireNamespace("sinkr", quietly = T))
      stop("Package 'sinkr' is required to run this function.\nYou can install from GitHub repository using the commands:\nlibrary(devtools)\ninstall_github('marchtaylor/sinkr')")
    
    # perform gap-filling
    if(verbose){
      message("Perform gap-filling")
    }
    rast_gapfilled <- sinkr::dineof(Xo=as.matrix(getValues(x@raster)), n.max=20, method="svds")
    # replace raster values with the gapfilled values
    values(x@raster) <- rast_gapfilled$Xa
    # mask result using input mask
    masked_pixels <- as.vector(as.numeric(which(seq(1:ncell(x@raster)) %in% na_index_mask == FALSE)))
    x@raster[masked_pixels] <- NaN
  }
  
  #################
  # Perform gap-filling using 'na.interpolation' method from package 'imputeTS'
  if(method %in% c("linear", "spline", "stine")){
    if(!requireNamespace("imputeTS", quietly = TRUE))
      stop("Package 'imputeTS' is required to run this function.\nYou can install using the command:\ninstall.packages('imputeTS', dependencies=TRUE)")
    
    # convert input raster to matrix
    rast_matrix <- as.matrix(x@raster[na_index_mask])
    
    # define interpolation function using selected method
    if(method %in% c("linear")){
      na.interpolation.par <- function(rmatrice){
        output <- imputeTS::na.interpolation(rmatrice, option="linear")
        return(output)
      }
    }
    
    if(method %in% c("spline")){
      na.interpolation.par <- function(rmatrice){
        output <- imputeTS::na.interpolation(rmatrice, option="spline")
        return(output)
      }
    }
    
    if(method %in% c("stine")){
      na.interpolation.par <- function(rmatrice){
        output <- imputeTS::na.interpolation(rmatrice, option="stine")
        return(output)
      }
    }
    
    # ### this part can be shortened if the line 'parallel::clusterExport(cl, varlist=c("method"), envir=environment())' works properly
    # na.interpolation.par <- function(rmatrice, method){
    #   output <- imputeTS::na.interpolation(rmatrice, option=method)
    #   return(output)
    # }
    
    # Perform gap-filling
    if(verbose){
      message("Perform gap-filling")
    }
    if(cores>1){
      if(!requireNamespace("parallel", quietly = TRUE)){
        warning("Package 'parallel' is required to run this function.\nGoing on using one core")
        matrix_gapfilled <- apply(X=rast_matrix, MARGIN=1, FUN=na.interpolation.par)
      } else {
        # check that machine has available cores
        if(cores > as.integer(parallel::detectCores())){
          cores <- as.integer(parallel::detectCores())
        }
        cl <- parallel::makePSOCKcluster(cores)
        #parallel::clusterExport(cl, varlist=c("method"), envir=environment())
        matrix_gapfilled <- parallel::parApply(cl=cl, X=rast_matrix, MARGIN=1, FUN=na.interpolation.par)
        parallel::stopCluster(cl)
      }
    } else {
      matrix_gapfilled <- apply(X=rast_matrix, MARGIN=1, FUN=na.interpolation.par)
    }
    
    # assemble results back to raster
    x@raster[na_index_mask] <- t(matrix_gapfilled)
    
    ### alternative method
    # rast_matrix_full <- matrix(data=NA, nrow=as.integer(ncell(x)), ncol=nlayers(x))
    # rast_matrix_full[na_index_mask,] <- t(matrix_gapfilled)
    # rast_stack <- stack(x@raster)
    # values(rast_stack) <- as.matrix(rast_matrix_full)
    # values(x@raster) <- getValues(rast_stack)
  }
  
  # return result
  return(x)
  
  # stop cluster
  on.exit(stopCluster(cl))
}
