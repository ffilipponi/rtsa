# title         : Mann-Kendall trend test of raster time series
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
# 
#' @title Raster time series Mann-Kendall trend test
#' 
#' @description This function conducts Mann-Kendall trend test from raster time series
#' using "Kendall" package especially designed to handle gappy time series.
#' 
#' @param rasterts Input raster time series as \code{\linkS4class{RasterStackTS}} or \code{\linkS4class{RasterBrickTS}} object.
#' @param rastermask Either a \code{\linkS4class{RasterLayer}} or "compute". Raster layer to use as a mask. When "compute" 
#' is set raster mask is computed to remove all pixels with incomplete time series.
#' @param gapfill Character. Defines the algorithm to be used to interpolate pixels with incomplete temporal profiles. 
#' Accepts argument supported as method in function \code{\link[rtsa]{rtsa.gapfill}}.
#' @param cores Integer. Defines the number of CPU to be used for multicore processing. Default to "1" core for 
#' singlecore processing.
#' @param ... Additional arguments to be passed through to function \code{\link[Kendall]{MannKendall}}.
#' 
#' @return Object of class \code{\link{MKstack-class}} containing the following components:
#' \tabular{rll}{
#' \tab \code{tau} \tab Kendall tau statistic\cr
#' \tab \code{pvalue} \tab Kendall two-sided p-value\cr
#' \tab \code{score} \tab Kendall Score\cr
#' \tab \code{variance} \tab Variance of Kendall Score\cr
#' }
#' 
#' @details 
#' 
#' @author Federico Filipponi
#' 
#' @references
#' 
#' Mann, H.B. (1945). Non-parametric tests against trend. Econometrica, 13, 163-171.
#' Kendall, M.G. (1975). Rank Correlation Methods, 4th edition. Charles Griffin, London.
#' Gilbert, R.O. (1987) . Statistical Methods for Environmental Pollution Monitoring. Wiley, NY.
#' Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.
#' Hipel, K.W. and McLeod, A.I., (2005). Time Series Modelling of Water Resources and Environmental Systems. 
#' Electronic reprint of our book orginally published in 1994. \href{http://www.stats.uwo.ca/faculty/aim/1994Book/}{book}
#' 
#' @keywords Mann-Kendall trend test time series analysis
#' 
#' @seealso \code{\link[Kendall]{MannKendall}}, \code{\link[rtsa]{rtsa.stl}}, \code{\link[rtsa]{rtsa.seas}}, \code{\link[rtsa]{rtsa.gapfill}}
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
#' # compute Mann-Kendall trend test
#' MannKendall_result <- rtsa.mk(rasterts=rasterts, rastermask=raster_mask)
#' # compute Mann-Kendall trend test using multiple cores on monthly time series
#' ### create monthly averages
#' rasterts_monthly_mean <- apply.monthly(rasterts, mean)
#' MannKendall_monhtly_result <- rtsa.mk(rasterts=rasterts_monthly_mean, rastermask=raster_mask)
#' }
#' 
#' @import raster
#' @import rts
#' @import parallel
#' @import doParallel
#' @importFrom Kendall MannKendall
#' @importFrom Kendall SeasonalMannKendall
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
#' @export

# define rsta.stl function
rtsa.mk <- function(rasterts, rastermask=NULL, gapfill="none", cores=1L, verbose=FALSE){
  
  # set environment
  options(scipen=7)
  
  # check if gapfill argument is set
  if(!(gapfill %in% c("none", "dineof", "linear", "spline", "stine")))
    stop("'gapfill' argument must be one of the following options: 'none', 'dineof', 'linear', 'spline', 'stine'")
  
  # require the 'stlplus' package
  if(!requireNamespace("Kendall", quietly = TRUE))
    stop("Package 'Kendall' is required to run this function.\nYou can install using the command:\ninstall.packages('Kendall', dependencies=TRUE)")

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
  
  # set 'mk_seasonal' variable based on periodicity
  if(seasonal_periodicity == 12){
    mk_seasonal <- TRUE
  } else {
    mk_seasonal <- FALSE
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
      message("Raster time series does not contain NA values after masking: OK.\nGap-filling will be not performed before the Mann-Kendall trend test computation") 
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
    }
  }
  
  # check number of NAs in time series

  if(na_check){
    # define 'countValid' function
    countValid <- function(x){
      output <- as.integer(length(which(!is.na(x))))
      return(output)
    }
    
    # count consecutive NAs
    if(cores>1){
      if(!requireNamespace("doParallel", quietly = TRUE)){
        warning("Package 'doParallel' is required to run this function.\nYou can install using the command:\ninstall.packages('doParallel', dependencies=TRUE)\nGoing on using one core")
        res_countValid <- apply(X=matrice, MARGIN=2, FUN=countValid)
      } else {
        # check that machine has available cores
        if(cores > as.integer(parallel::detectCores())){
          cores <- as.integer(parallel::detectCores())
        }
        # start cluster for parallel computation
        cl <- parallel::makePSOCKcluster(cores)
        doParallel::registerDoParallel(cl)
        # perform count of valid pixels
        res_countValid <- parallel::parCapply(cl=cl, x=matrice, FUN=countValid)
        # stop cluster
        parallel::stopCluster(cl)
      }
    } else {
      res_countValid <- apply(X=matrice, MARGIN=2, FUN=countValid)
    }
    
    # stop cluster
    #parallel::stopCluster(cl) ### check why it is not stopping inside the function
    
    # check if there are pixel profiles exceeding minimum number of valid observations for function 'MannKendall' (5)
    min_res_countValid <- min(res_countValid)
    if(min_res_countValid < 5){
      na_index_mask <- na_index_mask[-(which(min_res_countValid < 5))]
      if(length(na_index_mask) == 0){
        stop("Raster time series contain too many missing values and no pixel has a valid temporal profile.\nConsider the use of argument 'gapfill' to perform gapfilling before STL computation")
      }
      matrice <- as.matrix(rasterts[na_index_mask])
      if(verbose){
        message("Raster mask has been refined to deal with minumum number of time series oservations (minimum: 5 / found: ", min_res_countValid, ")")
        message("Final masked raster time series has ", ncol(matrice), " pixels and ", nrow(matrice), " temporal observations")
      }
    }
    
  }
  
  # compute Mann-Kendall test using 'Kendall' package

  # define parallel 'rtsa.mkpar' function
  # create time variable for the computation
  rtime <- as.integer(time(rasterts@time))
  
  # define Mann-Kendall function for parallel computing
  if(mk_seasonal){
    rtsa.mkpar <- function(x){
      tprofile <- xts::xts(as.vector(x), as.Date(rtime, origin="1970-01-01", tz="GMT"))
      tprofile2 <- stats::ts(data=as.vector(tprofile), start=c(as.integer(format(time(tprofile[1]), "%Y")), as.integer(format(time(tprofile[1]), "%m")), as.integer(format(time(tprofile[1]), "%d"))), end=c(as.integer(format(time(tprofile[length(tprofile)]), "%Y")), as.integer(format(time(tprofile[length(tprofile)]), "%m")), as.integer(format(time(tprofile[length(tprofile)]), "%d"))), deltat=deltats)
      mk_result <- Kendall::SeasonalMannKendall(tprofile2)
      output <- as.vector(c(as.double(mk_result$tau), as.double(mk_result$sl), as.integer(mk_result$S), as.integer(mk_result$varS)))
      return(output)
    }
  } else {
    rtsa.mkpar <- function(x){
      tprofile <- xts::xts(as.vector(x), as.Date(rtime, origin="1970-01-01", tz="GMT"))
      tprofile2 <- ts(data=as.vector(tprofile), start=c(as.integer(format(time(tprofile[1]), "%Y")), as.integer(format(time(tprofile[1]), "%m")), as.integer(format(time(tprofile[1]), "%d"))), end=c(as.integer(format(time(tprofile[length(tprofile)]), "%Y")), as.integer(format(time(tprofile[length(tprofile)]), "%m")), as.integer(format(time(tprofile[length(tprofile)]), "%d"))), deltat=deltats)
      mk_result <- Kendall::MannKendall(tprofile2)
      output <- as.vector(c(as.double(mk_result$tau), as.double(mk_result$sl), as.integer(mk_result$S), as.integer(mk_result$varS)))
      return(output)
    }
  }

  if(verbose){
    message(paste(c("Mann-Kendall trend test started at: "), Sys.time(), sep=""))
  }
  ptm <- proc.time()
  
  if(cores>1){
    if(!requireNamespace("doParallel", quietly = TRUE)){
      warning("Package 'doParallel' is required to run this function.\nYou can install using the command:\ninstall.packages('doParallel', dependencies=TRUE)\nGoing on using one core")
      res_mk <- apply(X=matrice, MARGIN=2, FUN=rtsa.mkpar)
    } else {
      # check that machine has available cores
      if(cores > as.integer(parallel::detectCores())){
        cores <- as.integer(parallel::detectCores())
      }
      if(verbose){
        message("Running Mann-Kendall trend test analysis using ", cores, " cores")
      }
      # run 'stl' computation in parallel
      # start cluster for parallel computation
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)
      # export 'stlplus' arguments to the cluster environment
      parallel::clusterExport(cl=cl, varlist=c("rtime", "deltats"), envir=environment())
      
      # perform 'stl' computation
      res_mk <- parallel::parCapply(cl=cl, x=matrice, FUN=rtsa.mkpar)
      
      # stop cluster
      parallel::stopCluster(cl)
      }
    } else {
      
      # run function using one core
      res_mk <- apply(X=matrice, MARGIN=2, FUN=rtsa.mkpar)
    }
  
  # stop cluster
  #parallel::stopCluster(cl) ### check why it is not stopping inside the function
  
  if(verbose){
    message(paste(c("Mann-Kendall trend test computation ended at: "), Sys.time(), sep=""))
    message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))
  }
  
  # create output object
  # assemble results to matrix object
  res_mk_matrix <- matrix(data=res_mk, nrow=ncol(matrice), ncol=4, byrow=TRUE)
  # create empty 'rts' object to store results
  void_raster <- raster(rasterts@raster[[1]])
  values(void_raster) <- NA
  
  # assemble results in a object of class 'STDstack' (Seasonal Trend Decomposition)
  mkreturn <- new("MKstack")

  void_raster[na_index_mask] <- as.vector(res_mk_matrix[,1])
  mkreturn@tau <- void_raster
  values(void_raster) <- NA
  
  void_raster[na_index_mask] <- as.vector(res_mk_matrix[,2])
  mkreturn@pvalue <- void_raster
  values(void_raster) <- NA
  
  void_raster[na_index_mask] <- as.vector(res_mk_matrix[,3])
  mkreturn@score <- void_raster
  values(void_raster) <- NA
  
  void_raster[na_index_mask] <- as.vector(res_mk_matrix[,4])
  mkreturn@variance <- void_raster

  # return function result
  return(mkreturn)
  
  # stop cluster
  on.exit(stopCluster(cl))
}
