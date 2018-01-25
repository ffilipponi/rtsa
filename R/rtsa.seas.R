# title         : Seasonal Adjustment analysis of raster time series
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
# 
#' @title Raster time series Seasonal Adjustment analysis using X-11 and X-13ARIMA-SEATS methods
#' 
#' @description This function conducts a Seasonal Adjustment analysis using X-11 and X-13ARIMA-SEATS methods 
#' from monthly raster time series using "seasonal" and "x13binary" packages.
#' 
#' @param rasterts Input raster time series as \code{\linkS4class{RasterStackTS}} or \code{\linkS4class{RasterBrickTS}} object.
#' @param rastermask Either a \code{\linkS4class{RasterLayer}} or "compute". Raster layer to use as a mask. When "compute" 
#' is set raster mask is computed to remove all pixels with incomplete time series.
#' @param method Character. Defines the method to be used for the Seasonal Adjustment analysis. Accepts argument 'x11' (default) or 'x13'.
#' @param gapfill Character. Defines the algorithm to be used to interpolate pixels with incomplete temporal profiles. 
#' Accepts argument supported as method in function \code{\link[rtsa]{rtsa.gapfill}}.
#' @param cores Integer. Defines the number of CPU to be used for multicore processing. Default to "1" core for 
#' singlecore processing.
#' @param only.statistics Logical. If TRUE returns only the statistics from seasonal, trend and remainder components.
#' @param keep.original Logical. If TRUE returns the original raster time series values in the 'rts' slot of \code{\link{STDstack-class}} object.
#' @param ... Additional arguments to be passed through to function \code{\link[seasonal]{seas}}.
#' 
#' @return Object of class \code{\link{STDstack-class}} containing the following components:
#' \tabular{rll}{
#' \tab \code{std} \tab Seasonal Trend Decomposition method used\cr
#' \tab \code{mask} \tab Final raster mask of computed pixels as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{seasonal_amplitude} \tab Amplitude of seasonal component (statistic) as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{seasonal_amplitude_stdev} \tab Standard deviation computed from the amplitude of seasonal component (statistic) as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{trend_slope} \tab Trend slope computed from trend component (yearly statistic) as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{residual_stdev} \tab Standard deviation computed from the remainder component (statistics) as \code{\linkS4class{RasterLayer}} object\cr
#' \tab \code{rts} \tab Input raster time series as \code{\linkS4class{RasterBrickTS}} object (only returned if \code{keep.original = TRUE})\cr
#' \tab \code{seasonality} \tab Seasonal component as \code{\linkS4class{RasterBrickTS}} object\cr
#' \tab \code{trend} \tab Trend component as \code{\linkS4class{RasterBrickTS}} object\cr
#' \tab \code{seasonaladjtrend} \tab Seasonal adjusted trend component as \code{\linkS4class{RasterBrickTS}} object\cr
#' \tab \code{remainder} \tab Remainder component as \code{\linkS4class{RasterBrickTS}} object\cr
#' }
#' 
#' @details 
#' 
#' @author Federico Filipponi
#' 
#' @references
#' 
#' Shiskin, J., Young, A.H., Musgrave, J.C. (1967). The X-11 variant of the Census Method II 
#' seasonal adjustment program. Technical Paper No. 15, U.S. Department of Commerce, U. S. Census Bureau.
#' 
#' Dagum, E.B. (1978). Modelling, forecasting and seasonally adjusting economic time series with the X-11 ARIMA method. 
#' Journal of the Royal Statistical Society. Series D (The Statistician), 27(3/4), 203-216.
#' 
#' Comprehensive list of R examples from the X-13ARIMA-SEATS manual: 
#' \href{http://www.seasonal.website/examples.html}{manual}
#' 
#' Official X-13ARIMA-SEATS manual: 
#' \href{https://www.census.gov/ts/x13as/docX13ASHTML.pdf}{pdf}
#' 
#' @keywords x-11 x-13ARIMA-SEATS seasonal adjustment time series analysis
#' 
#' @seealso \code{\link[seasonal]{seas}}, \code{\link[stlplus]{stlplus}}, \code{\link[rtsa]{rtsa.seas}}, \code{\link[rtsa]{rtsa.gapfill}}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}}
#' 
#' @examples
#' \dontrun{
#' ## create raster time series using the 'pacificSST' data from 'remote' package
#' require(remote)
#' 
#' data(pacificSST)
#' pacificSST[which(getValues(pacificSST == 0))] <- NA # set NA values
#' # subset input for faster processing
#' pacificSST_clip <- crop(pacificSST, extent(260, 290, -15, 15))
#' # create rts object
#' rasterts <- rts(pacificSST_clip, seq(as.Date('1982-01-15'), as.Date('2010-12-15'), 'months'))
#' 
#' ## generate raster mask
#' raster_mask <- pacificSST_clip[[1]] # create raster mask
#' names(raster_mask) <- "mask"
#' values(raster_mask) <- 1 # set raster mask values
#' raster_mask[which(is.na(getValues(pacificSST_clip[[1]])))] <- 0 # set raster mask values
#'
#' ## compute Seasonal Adjustment analysis
#' # use 'x11' (X-11) method to compute only adjusted seasonal trend decomposition statistics
#' std_x11_stats <- rtsa.seas(rasterts, rastermask=raster_mask, method="x11", only.statistics=TRUE)
#' # use 'x11' (X-11) method
#' std_x11_result <- rtsa.seas(rasterts=rasterts, rastermask=raster_mask, method="x11")
#' # use 'x13' (X-13ARIMA-SEATS) method
#' std_x13_result <- rtsa.seas(rasterts=rasterts, rastermask=raster_mask, method="x13")
#' # use 'x11' (X-11) method with multiple cores support and returning the original raster values
#' std_x11_res <- rtsa.seas(rasterts=rasterts, rastermask=raster_mask, cores=2, keep.original=TRUE)
#' }
#' 
#' @import raster
#' @import rts
#' @import parallel
#' @import doParallel
#' @import x13binary
#' @importFrom seasonal seas
#' @importFrom xts xts
#' @importFrom xts periodicity
#' @importFrom xts apply.yearly
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
rtsa.seas <- function(rasterts, rastermask=NULL, method="x11", gapfill="none", cores=1L, only.statistics=FALSE, keep.original=FALSE, verbose=FALSE){
  
  if(!(method %in% c("x11", "x13")))
    stop("'method' argument must be one of the following options: 'x11', 'x13'")
  
  # check if gapfill argument is set
  if(!(gapfill %in% c("none", "dineof", "gapfill", "linear", "spline", "stine")))
    stop("'gapfill' argument must be one of the following options: 'none', 'dineof', 'gapfill', 'linear', 'spline', 'stine'")
  
  # require the 'stlplus' package
  if(!requireNamespace("seasonal", quietly = TRUE))
    stop("Package 'seasonal' is required to run this function.\nYou can install using the command:\ninstall.packages('seasonal', dependencies=TRUE)")

  if(method == "x13"){
    # require the 'x13binary' package
    if(!requireNamespace("x13binary", quietly = TRUE))
      stop("Package 'x13binary' is required to run this function.\nYou can install using the command:\ninstall.packages('x13binary', dependencies=TRUE)")
  }
  
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
  
  # check length of time series and periodicity
  if(length(rasterts@time) > 780)
    stop("Input raster time series should have a maximum length of 780 observations.\nActual length is ", length(rasterts@time), " observations")
  
  if(seasonal_periodicity != 12) ### check if also > 12 would be feasible for seasonal time series
    stop("Periodicity of raster time series should be monthly")

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
      warning("Raster time series still contain NA values after masking.\nChange input raster mask or consider the use of the available options\nfor 'gapfill' argument to select a gap-filling method before the seasonal adjustment computation.\nSeasonal adjustment computation does not support gappy raster time series")
    }
  } else {
    if(verbose){
      message("Raster time series does not contain NA values after masking: OK.\nGap-filling will be not performed before the seasonal adjustment analysis computation") 
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
  }

  # compute Seasonal Trend Decomposition with seasonal adjustment using 'x11' or 'x13-arima-seats' using 'seasonal' and 'x13binary' packages'

  # define parallel 'rtsa.seaspar' function
  
  # create time variable for the computation
  rtime <- as.integer(time(rasterts@time))
  
  # define seas function based on method
  if(method == "x11"){
    rtsa.seasfun <- function(x){
      result <- tryCatch(expr=seasonal::seas(x, x11=""), error=function(e) NA)
    }
  }
  
  if(method == "x13"){
    rtsa.seasfun <- function(x){
      result <- tryCatch(expr=seasonal::seas(x), error=function(e) NA)
    }
  }

  rtsa.seaspar <- function(x){
    # create temporal object
    tprofile <- xts::xts(as.vector(x), as.Date(rtime, origin="1970-01-01", tz="GMT"))

    # convert xts to ts object (required by seasonal)
    tprofile2 <- stats::ts(data=as.vector(tprofile), start=c(as.integer(format(time(tprofile[1]), "%Y")), as.integer(format(time(tprofile[1]), "%m")), as.integer(format(time(tprofile[1]), "%d"))), end=c(as.integer(format(time(tprofile[length(tprofile)]), "%Y")), as.integer(format(time(tprofile[length(tprofile)]), "%m")), as.integer(format(time(tprofile[length(tprofile)]), "%d"))), deltat=deltats) 

    # compute seasonal trend adjustment analysis
    res_seas <- suppressMessages(rtsa.seasfun(tprofile2))
    
    # manage errors generated in function 'seas'
    if(is.na(res_seas)[1]){
      if(only.statistics){
        output <- as.vector(rep(NA, as.integer(4)))
      } else {
        output <- as.vector(rep(NA, as.integer(length(x)*4+4)))
      }
    } else {
      if(length(res_seas$data[1,])<5){
        if(only.statistics){
          output <- as.vector(rep(NA, as.integer(4)))
        } else {
          output <- as.vector(rep(NA, as.integer(length(x)*4+4)))
        }
      } else {
        
        # compute statistics from stlplus results
        tseasonal <- xts::xts(as.vector(res_seas$data[,2]), as.Date(rtime, origin="1970-01-01", tz="GMT"))
        seasonality_amplitude <- as.double(mean(xts::apply.yearly(tseasonal, max) - xts::apply.yearly(tseasonal, min)))
        seasonality_amplitude_stdev <- as.double(stats::sd(xts::apply.yearly(tseasonal, max) - xts::apply.yearly(tseasonal, min)))
        
        Xmin_pos <- 1
        Xmax_pos <- length(tprofile2)
        Xmin <- time(tprofile2)[Xmin_pos]
        Xmax <- time(tprofile2)[Xmax_pos]
        Ymin <- res_seas$data[,4][Xmin_pos]
        Ymax <- res_seas$data[,4][Xmax_pos]
        
        trend_slope <- as.double((Ymax-Ymin)/(Xmax-Xmin))
        
        remainder_stdev <- as.double(sd(res_seas$data[,5], na.rm=TRUE))
        
        # create function output
        if(only.statistics){
          output <- as.vector(c(seasonality_amplitude, seasonality_amplitude_stdev, trend_slope, remainder_stdev))
        } else {
          output <- as.vector(c(seasonality_amplitude, seasonality_amplitude_stdev, trend_slope, remainder_stdev, res_seas$data[,2], res_seas$data[,3], res_seas$data[,4], res_seas$data[,5]))
        }
      }
    }
    return(output)
  }

  if(verbose){
    message(paste(paste("Seasonal adjustment analysis using method '", method, "' started at: ", sep=""), Sys.time(), sep=""))
    message("It may take a long time ...")
  }
  ptm <- proc.time()
  
  if(cores>1){
    if(!requireNamespace("doParallel", quietly = TRUE)){
      warning("Package 'doParallel' is required to run this function.\nYou can install using the command:\ninstall.packages('doParallel', dependencies=TRUE)\nGoing on using one core")
      res_seas <- apply(X=matrice, MARGIN=2, FUN=rtsa.seaspar)
    } else {
      # check that machine has available cores
      if(cores > as.integer(parallel::detectCores())){
        cores <- as.integer(parallel::detectCores())
      }
      if(verbose){
        message("Running Seasonal Adjustment analysis using ", cores, " cores")
      }
      # run 'seasonal adjustment' computation in parallel
      # start cluster for parallel computation
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)
      # export 'seasonal' arguments to the cluster environment
      parallel::clusterExport(cl=cl, varlist=c("rtime", "deltats", "rtsa.seasfun", "only.statistics"), envir=environment())
      
      # perform 'seasonal adjustment' computation
      res_seas <- parallel::parCapply(cl=cl, x=matrice, FUN=rtsa.seaspar)
      
      # stop cluster
      parallel::stopCluster(cl)
      }
    } else {
      
      # run function using one core
      res_seas <- apply(X=matrice, MARGIN=2, FUN=rtsa.seaspar)
    }
  
  # stop cluster
  #parallel::stopCluster(cl) ### check why it is not stopping inside the function
  
  if(verbose){
    message(paste(c("Seasonal adjustment computation ended at: "), Sys.time(), sep=""))
    message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))
  }
  
  # create output object
  # assemble results to matrix object
  if(only.statistics){
    res_seas_matrix <- matrix(data=res_seas, nrow=ncol(matrice), ncol=4L, byrow=TRUE)
  } else {
    res_seas_matrix <- matrix(data=res_seas, nrow=ncol(matrice), ncol=as.integer(4+(nrow(matrice)*4)), byrow=TRUE)
  }
  
  # assemble results in a object of class 'STDstack' (Seasonal Trend Decomposition)
  seasreturn <- new("STDstack")
  
  seasreturn@std <- as.character(method)
  
  # fill in with mask and statistics
  void_raster <- raster(rasterts@raster[[1]])
  values(void_raster) <- NA
  
  names(void_raster) <- validNames("seasonal_amplitude")
  void_raster[na_index_mask] <- as.vector(res_seas_matrix[,1])
  seasreturn@seasonal_amplitude <- void_raster
  values(void_raster) <- NA
  
  names(void_raster) <- validNames("seasonal_amplitude_stdev")
  void_raster[na_index_mask] <- as.vector(res_seas_matrix[,2])
  seasreturn@seasonal_amplitude_stdev <- void_raster
  values(void_raster) <- NA
  
  names(void_raster) <- validNames("trend_slope")
  void_raster[na_index_mask] <- as.vector(res_seas_matrix[,3])
  seasreturn@trend_slope <- void_raster
  values(void_raster) <- NA
  
  names(void_raster) <- validNames("remainder_stdev")
  void_raster[na_index_mask] <- as.vector(res_seas_matrix[,4])
  seasreturn@remainder_stdev <- void_raster

  # store final mask
  names(void_raster) <- validNames("mask")
  na_index_mask_final <- which(!is.na(getValues(seasreturn@trend_slope)))

  values(void_raster) <- 0
  void_raster[na_index_mask_final] <- 1
  seasreturn@mask <- void_raster
  
  if(!only.statistics){
    # create empty 'rts' object to store results
    void_rasterts <- rasterts
    values(void_rasterts@raster) <- NA
    #stl_seasonality <- void_rasterts
    void_rasterts@raster[na_index_mask] <- as.vector(res_seas_matrix[,5:(nrow(matrice)+4)])
    seasreturn@seasonality <- void_rasterts
    values(void_rasterts@raster) <- NA
    
    void_rasterts@raster[na_index_mask] <- as.vector(res_seas_matrix[,((nrow(matrice)+5):((nrow(matrice)*2)+4))])
    seasreturn@trend <- void_rasterts
    values(void_rasterts@raster) <- NA
    
    void_rasterts@raster[na_index_mask] <- as.vector(res_seas_matrix[,(((nrow(matrice)*2)+5):((nrow(matrice)*3)+4))])
    seasreturn@seasonaladjtrend <- void_rasterts
    values(void_rasterts@raster) <- NA
    
    void_rasterts@raster[na_index_mask] <- as.vector(res_seas_matrix[,(((nrow(matrice)*3)+5):((nrow(matrice)*4)+4))])
    seasreturn@remainder <- void_rasterts
    }

  if(keep.original){
    seasreturn@rts <- rasterts
  }

  # return function result
  return(seasreturn)
  
  # stop cluster
  on.exit(stopCluster(cl))
}
