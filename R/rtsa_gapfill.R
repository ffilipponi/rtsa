# title         : Gap-filling raster time series
# Date          : Sep 2017
# Version       : 0.1
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
#
#' @title Raster time series gap-filling
#' 
#' @description This function perform gap-filling of gappy raster time series
#' 
#' @param x Input raster time series as \code{\link{RasterStackTS}} or \code{\link{RasterBrickTS}} object
#' @param rastermask A \code{\link{RasterLayer}} to use as a mask. If not set
#' a raster mask is computed to remove all pixels with less than two values in temporal profiles
#' @param method Character. Defines the algorithm to be used to interpolate pixels with incomplete temporal profiles. 
#' Accepts the following input:
#' \tab \code{"linear"} \tab for linear interpolation in \code{\link[imputeTS]{na.interpolation}} using \code{\link[stats]{approxfun}}
#' \tab \code{"spline"} \tab for spline interpolation in \code{\link[imputeTS]{na.interpolation}} using \code{\link[stats]{splinefun}}
#' \tab \code{"stine"} \tab for stine interpolation in \code{\link[imputeTS]{na.interpolation}} using \code{\link[stinepack]{stinterp}}
#' \tab \code{"dineof"} \tab for dineof interpolation using \code{\link[sinkr]{dineof}}
#' \tab \code{"gapfill"} \tab for gapfill interpolation using \code{\link[gapfill]{gapfill}}
#' @param codes Integer. Defines the number of CPU to be used for multicore processing. Default to "1" core for 
#' singlecore processing.
#' @param ... Additional arguments
#' 
#' @return Object of class \code{\link{RasterBrickTS}} with gap-filled pixels
#' 
#' @details 
#' 
#' @author Federico Filipponi
#' 
#' @references
#' 
#' @keywords time series analysis gap-filling
#' 
#' @seealso \code{\link[imputeTS]{na.interpolation}} \code{\link[sinkr]{dineof}}, \code{\link[gapfill]{gapfill}}, \code{\link[stats]{approxfun}}, \code{\link[stats]{splinefun}}, \code{\link[stinepack]{stinterp}}
#' 
#' @examples
#' #' \dontrun{
#' ## create raster time series using the 'pacificSST' data from 'remote' package
#' require(remote)
#' 
#' data(pacificSST)
#' pacificSST[which(getValues(pacificSST == 0))] <- NA # set NA values
#' rasterts <- rts(pacificSST, seq(as.Date('1982-01-15'), as.Date('2010-12-15'), 'months')) # create rts object
#' 
#' ## generate raster mask
#' raster_mask <- pacificSST[[1]] # create raster mask
#' values(raster_mask) <- 1 # set raster mask values
#' raster_mask[which(is.na(getValues(pacificSST[[1]])))] <- 0 # set raster mask values
#' 
#' ## randomly remove values from cells in rts object
#' frac_gaps <- 0.5 # the fraction of data with NaNs
#' temporal_cells <- as.integer(ncell(rasterts) * nlayers(rasterts)) # number of total cells in the rts
#' na_cells <- sort(unique(sample.int(temporal_cells, (temporal_cells * frac_gaps)))) # define random position of cells to be set to NaN
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
#' 
#' # using gapfill interpolation and raster mask
#' rasterts_gapfill <- rtsa.gapfill(rasterts_gappy, rastermask=raster_mask, method="gapfill")
#' }
#' 
#' @import sinkr gapfill imputeTS doParallel
#' 
#' @export

rtsa.gapfill <- function(x, rastermask=NULL, method, cores=1L, verbose=FALSE){
  
  if(!(class(x) %in% c("RasterStackTS", "RasterBrickTS")))
    stop("'x' argument must be an object of class 'RasterStackTS', 'RasterBrickTS'")
  
  if(!(method %in% c("dineof", "gapfill", "linear", "spline", "stine")))
    stop("'method' argument must be one of the following options: 'dineof', 'gapfill', 'linear', 'spline', 'stine'")
  
  #################
  # Process raster mask
  if(!(is.null(rastermask))){
    if(!(class(rastermask) %in% c("RasterLayer"))){
      warning("'rastermask' argument is not an object of class 'raster'.\nNo raster mask will be used to clip the input raster time series")
      na_index_mask <- as.vector(1:as.integer(ncell(x)))
    } else {
      # compare mask raster extent with the input raster time series
      mask_correspondence <- as.logical(compareRaster(x@raster, rastermask, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, rotation=TRUE, values=FALSE, stopiffalse=FALSE))
      #mask_correspondence <- TRUE ### for debug
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
      message("Mask raster not supplied. Masking pixel with less than two temporal observation")
    }
    #na_index_mask <- as.vector(1:as.integer(ncell(x)))
    
    # compute mask using pixel that have less than two values in temporal profiles
    # define function to mask incomplete pixel profiles
    generateMask <- function(x){
      m <- as.integer(length(which(is.na(x)))<=(length(x)-2))
      return(m)
    }
    # generate raster mask
    generated_mask <- as.vector(apply(X=as.matrix(getValues(x)), FUN=generateMask, MARGIN=1))
    ### this may be parallelized after setting a cluster using: parApply(cl=NULL, X=as.matrix(getValues(x)), FUN=generateMask, MARGIN=1)
    na_index_mask <- as.vector(which(generated_mask==1))
  }
  
  #################
  # Perform gap-filling using 'dineof' method from package 'sinkr'
  if(method == "dineof"){
    if(!require("sinkr")){
      stop("Package 'sinkr' is required to run this function.\nYou can install from GitHub repository using the commands:\nlibrary(devtools)\ninstall_github('marchtaylor/sinkr')")
    }
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
  # Perform gap-filling using 'gapfill' method from package 'gapfill'
  if(method == "gapfill"){
    if(!require("gapfill")){
      stop("Package 'gapfill' is required to run this function.\nYou can install using the command:\ninstall.packages('gapfill', dependencies=TRUE)")
    }
    ### add control requireing at least 4 images
    
    # extract year
    xts_year <- format(as.yearmon(index(x@time)), "%Y")
    # year range
    year_range <- as.integer((max(as.numeric(xts_year)) - min(as.numeric(xts_year)))+1)
    
    # continue only if number of observation years is greater than 1
    if(year_range <= 1){
      stop("Number of observation years should be greater than one.\nSelected another gapfilling method")
    } else {
      # year unique values
      year_unique <- unique(min(as.numeric(xts_year)):max(as.numeric(xts_year)))
      # extract doy
      xts_doy <- as.integer(format(as.POSIXct(index(x@time), origin="1970-01-01", tz="GMT"), "%j"))
      # extract doy range
      doy_range <- as.integer((max(as.numeric(xts_doy)) - min(as.numeric(xts_doy)))+1)
      # doy unique values
      doy_unique <- unique(min(as.numeric(xts_doy)):max(as.numeric(xts_doy)))
      
      # extract periodicity ### use periodicity to create
      # create data.frame for scale conversion
      #time_scale_conv <- data.frame(periodicity=c("yearly", "quarter", "monthly", "weeks", "days"), seq=c("years", "quarter", "months", "weeks", "days"))
      
      # extract time periodicity of raster time series
      time_periodicity <- periodicity(x@time)
      #time_scale <- time_periodicity$scale
      
      # create regularly spaced sequence of dates
      obs_time <- seq(time_periodicity$start, time_periodicity$end, by=time_periodicity$label)
      # clean 29 Feb if periodicity is different from days
      if(time_periodicity$label != "day"){
        obs_time_clean <- obs_time[format(obs_time, "%m %d") !=  "02 29"]
      }
      # extract doy range
      doy_range <- as.integer(as.integer(range(format(obs_time, "%j"))[2]) - as.integer(range(format(obs_time, "%j"))[1])+1)
      # clean 
      
      doy_unique <- sort(unique(as.integer(format(obs_time, "%j"))))
      # year range
      # year_range <- as.integer(as.integer(range(format(obs_time, "%Y"))[2]) - as.integer(range(format(obs_time, "%Y"))[1])+1)
      
      # create empty array to store raster time series
      arr <- array(dim=c(nrow(x), ncol(x), doy_range, year_range))
      # set array dimansion names
      dimnames(arr)[[1]] <- as.vector(sort(unique(coordinates(raster(x[[1]]))[,2])))
      dimnames(arr)[[2]] <- as.vector(sort(unique(coordinates(raster(x[[1]]))[,1])))
      #dimnames(arr)[[3]] <- as.vector(doy_unique) ### gives error
      dimnames(arr)[[4]] <- year_unique
      
      # fill in array with raster time series values
      for(y in min(as.numeric(xts_year)):max(as.numeric(xts_year))){
        for(d in min(as.numeric(xts_doy)):max(as.numeric(xts_doy))){
          layer_id <- as.integer(which(xts_year == y & xts_doy == d)[1])
          arr_y <- which(year_unique == y)
          arr_d <- which(doy_unique == d)
          if(!is.na(as.integer(which(xts_year == y & xts_doy == d)[1]))){
            arr[,,arr_d,arr_y] <- as.matrix(x@raster[[layer_id]])
          } else {
            arr[,,arr_d,arr_y] <- NaN
          }
        }
      }
      
      # get mask values
      na_index_matrix <- as.matrix(rastermask)
      # create arry for masking
      mask_array <- array(dim=c(x@raster@nrows, x@raster@ncols, doy_range, year_range))
      # set array dimansion names
      dimnames(mask_array)[[1]] <- as.vector(sort(unique(coordinates(raster(x[[1]]))[,2])))
      dimnames(mask_array)[[2]] <- as.vector(sort(unique(coordinates(raster(x[[1]]))[,1])))
      dimnames(mask_array)[[3]] <- as.vector(doy_unique)
      dimnames(mask_array)[[4]] <- year_unique
      # set values for the array
      for(y in min(as.numeric(xts_year)):max(as.numeric(xts_year))){
        for(d in min(as.numeric(xts_doy)):max(as.numeric(xts_doy))){
          #layer_id <- as.integer(which(xts_year == y & xts_doy == d)[1])
          arr_y <- which(year_unique == y)
          arr_d <- which(doy_unique == d)
          if(!is.na(as.integer(which(xts_year == y & xts_doy == d)[1]))){
            mask_array[,,arr_d,arr_y] <- as.logical(na_index_matrix)
          } else {
            mask_array[,,arr_d,arr_y] <- matrix(data=as.logical(0), nrow=x@raster@nrows, ncol=x@raster@ncols)
          }
        }
      }
      
      # Perform gap-filling
      if(verbose){
        message("Perform gap-filling. It may take a long time")
      }
      if(cores>1){
        if(!require(doParallel)){
          warning("Package 'doParallel' is required to run this function.\nYou can install using the command:\ninstall.packages('doParallel', dependencies=TRUE)\nGoing on using one core")
          arr_gapfilled <- gapfill::Gapfill(data=arr, subset=mask_array, dopar=FALSE)
        } else {
          # check that machine has available cores
          if(cores > as.integer(parallel::detectCores())){
            cores <- as.integer(parallel::detectCores())
          }
          cl <- makePSOCKcluster(cores)
          doParallel::registerDoParallel(cl)
          #cl <- parallel::makePSOCKcluster(cores)
          arr_gapfilled <- gapfill::Gapfill(data=arr, subset=mask_array, dopar=TRUE)
          stopCluster(cl)
        }
      } else {
        arr_gapfilled <- gapfill::Gapfill(data=arr, subset=mask_array, dopar=FALSE)
      }
      
      # assemble results back to raster
      for(y in min(as.numeric(xts_year)):max(as.numeric(xts_year))){
        for(d in min(as.numeric(xts_doy)):max(as.numeric(xts_doy))){
          layer_id <- as.integer(which(xts_year == y & xts_doy == d)[1])
          arr_y <- which(year_unique == y)
          arr_d <- which(doy_unique == d)
          if(!is.na(as.integer(which(xts_year == y & xts_doy == d)[1]))){
            x@raster[[layer_id]] <- t(as.matrix(arr_gapfilled$fill[,,arr_d,arr_y]))
          } else {
            #
          }
        }
      }
    }
  }
  
  #################
  # Perform gap-filling using 'na.interpolation' method from package 'imputeTS'
  if(method %in% c("linear", "spline", "stine")){
    if(!require("imputeTS")){
      stop("Package 'imputeTS' is required to run this function.\nYou can install using the command:\ninstall.packages('imputeTS', dependencies=TRUE)")
    }
    
    # # perform gap-filling directly on raster object
    # for(p in na_index_mask){
    #   x@raster[p] <- imputeTS::na.interpolation(x=as.vector(x@raster[p]), option=method)
    # }
    
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
      if(!require(parallel)){
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
}
