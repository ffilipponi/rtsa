# title         : Seasonal Trend Decomposition using Loess (STL) of raster time series
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
# 
#' @title Raster time series STL (Seasonal Trend Decomposition using Loess)
#' 
#' @description This function conducts a Seasonal Trend Decomposition using Loess (STL) from raster time series
#' using "stlplus" package especially designed to handle gappy time series.
#' 
#' @param rasterts Input raster time series as \code{\linkS4class{RasterStackTS}} or \code{\linkS4class{RasterBrickTS}} object.
#' @param rastermask Either a \code{\linkS4class{RasterLayer}} or "compute". Raster layer to use as a mask. When "compute" 
#' is set raster mask is computed to remove all pixels with incomplete time series.
#' @param gapfill Character. Defines the algorithm to be used to interpolate pixels with incomplete temporal profiles. 
#' Accepts argument supported as method in function \code{\link[rtsa]{rtsa.gapfill}}.
#' @param cores Integer. Defines the number of CPU to be used for multicore processing. Default to "1" core for 
#' singlecore processing.
#' @param n.p Integer. Argument to be passed to function \code{\link[stlplus]{stlplus}}. Periodicity of the seasonal component. 
#' Default to the frequency estimated from of the time series.
#' @param s.window Argument to be passed to function \code{\link[stlplus]{stlplus}}. Either the character string "periodic" (default) 
#' or the span (in lags) of the loess window for seasonal extraction, which should be odd.
#' @param t.window Integer. Argument to be passed to function \code{\link[stlplus]{stlplus}}. The span (in lags) of the the loess window for trend extraction, 
#' which should be odd. Default to the entire time series duration.
#' @param s.degree Integer. Argument to be passed to function \code{\link[stlplus]{stlplus}}. Degree of locally-fitted polynomial in trend extraction.
#' Should be 0, 1 (default) or 2.
#' @param t.degree Integer. Argument to be passed to function \code{\link[stlplus]{stlplus}}. Degree of locally-fitted polynomial in trend extraction.
#' Should be 0, 1 (default) or 2.
#' @param only.statistics Logical. If TRUE returns only the statistics from seasonal, trend and remainder components.
#' @param keep.original Logical. If TRUE returns the original raster time series values in the 'rts' slot of \code{\link{STDstack-class}} object.
#' @param ... Additional arguments to be passed through to function \code{\link[stlplus]{stlplus}}.
#' 
#' @return Object of class \code{\link{STDstack-class}} containing the following components:
#' \tabular{rll}{
#' \tab \code{std} \tab Seasonal Trend Decomposition method used\cr
#' \tab \code{mask} \tab Final raster mask of computed pixels as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{seasonal_amplitude} \tab Amplitude of seasonal component (statistic) as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{seasonal_amplitude_stdev} \tab Standard deviation computed from the amplitude of seasonal component (statistic) as \code{\linkS4class{RasterLayer}} object (only returned when running function \code{\link[rtsa]{rtsa.seas}})\cr
#' \tab \code{trend_slope} \tab Trend slope computed from trend component (yearly statistic) as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{residual_stdev} \tab Standard deviation computed from the remainder component (statistics) as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{rts} \tab Input raster time series as \code{\linkS4class{RasterBrickTS}} object (only returned if \code{keep.original = TRUE})\cr
#' \tab \code{seasonality} \tab Seasonal component as \code{\linkS4class{RasterBrickTS}} object\cr
#' \tab \code{trend} \tab Trend component as \code{\linkS4class{RasterBrickTS}} object\cr
#' \tab \code{seasonaladjtrend} \tab Seasonal adjusted trend component as \code{\linkS4class{RasterBrickTS}} object (only returned when running function \code{\link[rtsa]{rtsa.seas}})\cr
#' \tab \code{remainder} \tab Remainder component as \code{\linkS4class{RasterBrickTS}} object\cr
#' }
#' 
#' @details 
#' 
#' @author Federico Filipponi
#' 
#' @references
#' 
#' Cleveland, R.B., Cleveland, W.S., Terpenning, I. (1990). STL: A seasonal-trend decomposition procedure based on loess. 
#' Journal of Official Statistics, 6(1), 3.
#' 
#' Hafen, R.P. (2010). Local regression models: Advancements, applications, and new methods. 
#' West Lafayette, Indiana: Purdue University, PhD dissertation, pp. 279.
#' 
#' @keywords STL Loess Seasonal Trend decomposition time series analysis
#' 
#' @seealso \code{\link[stlplus]{stlplus}}, \code{\link[rtsa]{rtsa.seas}}, \code{\link[rtsa]{rtsa.gapfill}}, \code{\link[seasonal]{seas}}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}}
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
#' ## compute Seasonal Trend Decomposition analysis
#' # compute 'STL'
#' std_result <- rtsa.stl(rasterts=rasterts, rastermask=raster_mask)
#' # compute STL' with multiple cores support and returning the original raster values
#' std_result <- rtsa.stl(rasterts=rasterts, rastermask=raster_mask, cores=2, keep.original=TRUE)
#' # compute STL' using additional arguments to define loess window
#' std_result <- rtsa.stl(rasterts=rasterts, rastermask=raster_mask, n.p=12, t.window=48)
#' }
#' 
#' @import raster
#' @import rts
#' @import parallel
#' @import doParallel
#' @import stlplus
#' @importFrom stlplus stlplus
#' @importFrom xts xts
#' @importFrom xts periodicity
#' @importFrom parallel detectCores
#' @importFrom parallel makePSOCKcluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parCapply
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom methods new
#' @importFrom stats time
#' @importFrom stats ts
#' @importFrom stats sd
#' @export

# define rsta.stl function
rtsa.stl <- function(rasterts, rastermask=NULL, gapfill="none", cores=1L, n.p=NULL, s.window="periodic", t.window=NULL, s.degree=1L, t.degree=1L, only.statistics=FALSE, keep.original=FALSE, verbose=FALSE){
  
  # check if gapfill argument is set
  if(!(gapfill %in% c("none", "dineof", "linear", "spline", "stine")))
    stop("'gapfill' argument must be one of the following options: 'none', 'dineof', 'linear', 'spline', 'stine'")
  
  # require the 'stlplus' package
  if(!requireNamespace("stlplus", quietly = TRUE))
    stop("Package 'stlplus' is required to run this function.\nYou can install using the command:\ninstall.packages('stlplus', dependencies=TRUE)")

  # check if input file is an object of class 'RasterStackTS', 'RasterBrickTS'
  if(!(class(rasterts) %in% c("RasterStackTS", "RasterBrickTS")))
    stop("'rasterts' argument must be an object of class 'RasterStackTS' or 'RasterBrickTS'.\nUse 'rts()' function to generate 'rasterts' input")
  
  # set time periodicity
  # define table for periodicity to deltat conversion
  deltatable <- data.frame(scale=c("yearly", "quarterly", "monthly", "weekly", "daily", "hourly", "minute", "seconds"), 
                           deltat=c(1, 1.0/3, 1.0/12, 1.0/52.17857, 1.0/365, 1.0/8760, 1.0/525600, 1.0/31536000), 
                           periodicity=c(1, 3, 12, 52.17857, 365, 8760, 525600, 31536000))
  
  # convert periodicity to deltat
  periodicity_tprofile <- xts::periodicity(rasterts@time)$scale
  deltats <- as.double(deltatable$deltat[which(deltatable$scale == periodicity_tprofile)])
  seasonal_periodicity <- as.double(deltatable$periodicity[which(deltatable$scale == periodicity_tprofile)])

  # set NULL input arguments to defaults
  
  # define function to check input smoothing degree (from 'stlplus.R')
  rtsa.degcheck <- function(x){
    if(!all(x ==0 | x == 1 | x == 2)) stop("Smoothing degree must be 0, 1 or 2")
  }
  
  if(is.null(n.p)){
    n.p <- seasonal_periodicity
  }
  
  if(is.null(t.window)){
    t.window <- as.integer(length(rasterts@time))
  }
  
  if(t.window != as.integer(length(rasterts@time))){
    warning(paste("'t.window' size (", t.window, ") is not equal to the length of time series (", length(rasterts@time), ").\nThis may generate incorrect trend slope values", sep=""))
  }
  
  if(!(s.window == "periodic")){
    if(is.numeric(s.window)){
      if(!(as.logical(s.window %% 2))){
        stop("Argument 's.window' must be a odd number or 'periodic'")
      }
      if(s.window > t.window){
        stop("Argument 't.window' must be greater than 's.window'")
      }
    } else {
      stop("Argument 's.window' should be either numeric (odd) or 'periodic'")
    }
    # check smoothing degree values
    rtsa.degcheck(s.degree)
  }
  
  # check smoothing degree values
  rtsa.degcheck(t.degree)
  
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
  
  # check if there are NAs in the masked dataset
  na_check <- anyNA(matrice)
  if(na_check){
    if(gapfill == "none"){
      warning("Raster time series still contain NA values after masking.\nChange input raster mask or consider the use of the available options\nfor 'gapfill' argument to select a gap-filling method before the STL computation.\nGoing on however using gappy input raster time series for the STL computation")
    }
  } else {
    if(verbose){
      message("Raster time series does not contain NA values after masking: OK.\nGap-filling will be not performed before the STL computation") 
    }
  }
  
  # perform gap-filling
  if(na_check){
    if(gapfill != "none"){
      # check if gapfill argument is set
      if(!(gapfill %in% c("none", "dineof", "gapfill", "linear", "spline", "stine")))
        stop("'gapfill' argument must be one of the following options: 'none', 'dineof', 'gapfill', 'linear', 'spline', 'stine'")
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
      matrice <- as.matrix(rasterts[na_index_mask])
    }
  }
  
  # check number of consecutive NAs in time series
  
  if(na_check){
    # define 'consecutiveNAcount' function
    consecutiveNAcount <- function(x){
      lagsize <- as.integer(floor(n.p/2))
      conta <- rep(0, length(x))
      for(j in 1:length(x)){
        conta[j] <- sum(as.integer(is.na(x[max(1,(j-lagsize)):min(length(x),(j+lagsize))])))
        if(length(which(is.na(x))) == length(x)){
          conta[1] <- length(x)
        }
      }
      output <- max(as.vector(conta))
      return(output)
    }
    
    # count consecutive NAs
    if(cores>1){
      if(!requireNamespace("doParallel", quietly = TRUE)){
        warning("Package 'doParallel' is required to run this function.\nYou can install using the command:\ninstall.packages('doParallel', dependencies=TRUE)\nGoing on using one core")
        res_consNA <- apply(X=matrice, MARGIN=2, FUN=consecutiveNAcount)
      } else {
        # check that machine has available cores
        if(cores > as.integer(parallel::detectCores())){
          cores <- as.integer(parallel::detectCores())
        }
        # start cluster for parallel computation
        cl <- parallel::makePSOCKcluster(cores)
        doParallel::registerDoParallel(cl)
        # export 'n.p' lag size to the cluster environment
        parallel::clusterExport(cl=cl, "n.p", envir=environment())
        # perform 
        res_consNA <- parallel::parCapply(cl=cl, x=matrice, FUN=consecutiveNAcount)
        # stop cluster
        parallel::stopCluster(cl)
      }
    } else {
      res_consNA <- apply(X=matrice, MARGIN=2, FUN=consecutiveNAcount)
    }
    
    # stop cluster
    #parallel::stopCluster(cl) ### check why it is not stopping inside the function
    
    # check if there are pixel profile exceeding maximum number of allowed consecutive gaps
    max_res_consNA <- max(res_consNA)
    if(max_res_consNA >= n.p){
      na_index_mask <- na_index_mask[-(which(res_consNA >= n.p))]
      if(length(na_index_mask) == 0){
        stop("Raster time series contain too many missing values and no pixel has a valid temporal profile.\nConsider the use of argument 'gapfill' to perform gapfilling before STL computation")
      }
      matrice <- as.matrix(rasterts[na_index_mask])
      if(verbose){
        message("Raster mask has been refined to deal with maximum number of allowed consecutive gaps (allowed: ", n.p-1, " / found: ", max_res_consNA, ")")
        message("Final masked raster time series has ", ncol(matrice), " pixels and ", nrow(matrice), " temporal observations")
      }
    }
  }
  
  # compute Seasonal Trend Decomposition using Loess (STL) using 'stlplus' package

  # create time variable for the computation
  rtime <- as.integer(time(rasterts@time))
  
  # define parallel 'rtsa.stlpar' function
  rtsa.stlpar <- function(x){
    # create temporal object
    tprofile <- xts::xts(as.vector(x), as.Date(rtime, origin="1970-01-01", tz="GMT"))

    # convert xts to ts object (required by stlplus)
    tprofile2 <- stats::ts(data=as.vector(tprofile), start=c(as.integer(format(time(tprofile[1]), "%Y")), as.integer(format(time(tprofile[1]), "%m")), as.integer(format(time(tprofile[1]), "%d"))), end=c(as.integer(format(time(tprofile[length(tprofile)]), "%Y")), as.integer(format(time(tprofile[length(tprofile)]), "%m")), as.integer(format(time(tprofile[length(tprofile)]), "%d"))), deltat=deltats) 

    # compute STL
    stl_profile_unique <- stlplus::stlplus(as.vector(tprofile2), s.window=s.window, n.p=n.p, t.window=t.window, s.degree=s.degree, t.degree=t.degree)
    
    # compute statistics from stlplus results
    seasonality_amplitude <- diff(range(stl_profile_unique$data$seasonal))
    #seasonal_amplitude_value <- as.double(seasonality_amplitude + mean(x, na.rm=TRUE))
    
    Xmin_pos <- 1
    Xmax_pos <- length(tprofile2)
    Xmin <- time(tprofile2)[Xmin_pos]
    Xmax <- time(tprofile2)[Xmax_pos]
    Ymin <- stl_profile_unique$data$trend[Xmin_pos]
    Ymax <- stl_profile_unique$data$trend[Xmax_pos]
    
    trend_slope <- as.double((Ymax-Ymin)/(Xmax-Xmin))
    
    remainder_stdev <- stats::sd(stl_profile_unique$data$remainder, na.rm=TRUE)
    
    # create function output
    if(only.statistics){
      output <- as.vector(c(seasonality_amplitude, trend_slope, remainder_stdev))
    } else {
      output <- as.vector(c(seasonality_amplitude, trend_slope, remainder_stdev, stl_profile_unique$data$seasonal, stl_profile_unique$data$trend, stl_profile_unique$data$remainder))
    }
    
    return(output)
  }
  
  if(verbose){
    message(paste(c("STL analysis started at: "), Sys.time(), sep=""))
    message("STL parameter set: n.p=", n.p, ", s.window=", s.window, ", t.window=", t.window, ", s.degree=", s.degree, ", t.degree=", t.degree)
  }
  ptm <- proc.time()
  
  if(cores>1){
    if(!requireNamespace("doParallel", quietly = TRUE)){
      warning("Package 'doParallel' is required to run this function.\nYou can install using the command:\ninstall.packages('doParallel', dependencies=TRUE)\nGoing on using one core")
      res_stl <- apply(X=matrice, MARGIN=2, FUN=rtsa.stlpar)
    } else {
      # check that machine has available cores
      if(cores > as.integer(parallel::detectCores())){
        cores <- as.integer(parallel::detectCores())
      }
      if(verbose){
        message("Running STL analysis using ", cores, " cores")
      }
      # run 'stl' computation in parallel
      # start cluster for parallel computation
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)
      # export 'stlplus' arguments to the cluster environment
      parallel::clusterExport(cl=cl, varlist=c("rtime", "deltats", "n.p", "t.window", "s.window", "s.degree", "t.degree", "only.statistics"), envir=environment())
      
      # perform 'stl' computation
      res_stl <- parallel::parCapply(cl=cl, x=matrice, FUN=rtsa.stlpar)
      
      # stop cluster
      parallel::stopCluster(cl)
      }
    } else {
      
      # run function using one core
      res_stl <- apply(X=matrice, MARGIN=2, FUN=rtsa.stlpar)
    }
  
  # stop cluster
  #parallel::stopCluster(cl) ### check why it is not stopping inside the function
  
  if(verbose){
    message(paste(c("STL computation ended at: "), Sys.time(), sep=""))
    message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))
  }
  
  # create output object
  # assemble results to matrix object
  if(only.statistics){
    res_stl_matrix <- matrix(data=res_stl, nrow=ncol(matrice), ncol=3L, byrow=TRUE)
  } else {
    res_stl_matrix <- matrix(data=res_stl, nrow=ncol(matrice), ncol=as.integer(3+(nrow(matrice)*3)), byrow=TRUE)
  }
  
  # assemble results in a object of class 'STDstack' (Seasonal Trend Decomposition)
  stlreturn <- new("STDstack")
  
  stlreturn@std <- as.character("stlplus")
  
  # fill in with mask and statistics
  void_raster <- raster(rasterts@raster[[1]])
  
  names(void_raster) <- validNames("mask")
  values(void_raster) <- 0
  void_raster[na_index_mask] <- 1
  stlreturn@mask <- void_raster
  values(void_raster) <- NA
  
  names(void_raster) <- validNames("seasonal_amplitude")
  void_raster[na_index_mask] <- as.vector(res_stl_matrix[,1])
  stlreturn@seasonal_amplitude <- void_raster
  values(void_raster) <- NA
  
  names(void_raster) <- validNames("trend_slope")
  void_raster[na_index_mask] <- as.vector(res_stl_matrix[,2])
  stlreturn@trend_slope <- void_raster
  values(void_raster) <- NA
  
  names(void_raster) <- validNames("remainder_stdev")
  void_raster[na_index_mask] <- as.vector(res_stl_matrix[,3])
  stlreturn@remainder_stdev <- void_raster
  #values(void_raster) <- NA
  
  if(!only.statistics){
    # create empty 'rts' object to store results
    void_rasterts <- rasterts
    values(void_rasterts@raster) <- NA
    #stl_seasonality <- void_rasterts
    void_rasterts@raster[na_index_mask] <- as.vector(res_stl_matrix[,4:(nrow(matrice)+3)])
    stlreturn@seasonality <- void_rasterts
    values(void_rasterts@raster) <- NA
    
    void_rasterts@raster[na_index_mask] <- as.vector(res_stl_matrix[,((nrow(matrice)+4):((nrow(matrice)*2)+3))])
    stlreturn@trend <- void_rasterts
    values(void_rasterts@raster) <- NA
    
    void_rasterts@raster[na_index_mask] <- as.vector(res_stl_matrix[,(((nrow(matrice)*2)+4):((nrow(matrice)*3)+3))])
    stlreturn@remainder <- void_rasterts
  }
  
  if(keep.original){
    stlreturn@rts <- rasterts
  }

  # return function result
  return(stlreturn)
  
  # stop cluster
  on.exit(stopCluster(cl))
}
