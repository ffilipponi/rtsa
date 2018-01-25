# title         : Empirical Orthogonal Teleconnection analysis of raster time series
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
# 
#' @title EOT (Empirical Orthogonal Teleconnections analysis)
#' 
#' @description This function conducts an Empirical Orthogonal Teleconnection analysis (EOT) of a raster time series 
#' predictor using "remote" package especially designed to identify locations of enhanced potential and explain 
#' spatio-temporal variability within the same geographic field.
#' 
#' @param rasterts Input raster time series as \code{\linkS4class{RasterStackTS}} or \code{\linkS4class{RasterBrickTS}} object.
#' @param rastermask Either a \code{\linkS4class{RasterLayer}} or "compute". Raster layer to use as a mask. When "compute" 
#' is set raster mask is computed to remove all pixels with incomplete time series.
#' @param nu Numeric. Defines the number of EOTs to return. Defaults to return the first 2 EOT modes.
#' @param gapfill Character. Defines the algorithm to be used to interpolate pixels with incomplete temporal profiles
#' Accepts argument supported as method in function \code{\link[rtsa]{rtsa.gapfill}}.
#' @param predictor Character. Defines the predictor components to export from those available in Value section of \code{\link[remote]{eot}}.
#' Supports one or more of the following: 'all', 'r_predictor', 'rsq_predictor', 'rsq_sums_predictor', 'int_predictor', 'slp_predictor', 'p_predictor'.
#' @param ... Additional arguments to be passed through to function \code{\link[remote]{eot}}.
#' 
#' @return Object of class \code{\link{EOTstack-class}} containing the following components:
#' \tabular{rll}{
#' \tab \code{eot} \tab EOT temporal profiles corresponding to base point coordinates as \code{\linkS4class{xts}} object\cr
#' \tab \code{total_variance} \tab Numeric. Total explained variance of input raster time series by the entire set of computed EOTs\cr
#' \tab \code{explained_variance} \tab Numeric vector. Percentage of variance explained by each EOT mode with respect to the total variance of input raster time series\cr
#' \tab \code{r_predictor} \tab RasterBrick. Correlation coefficients between the base point and each pixel of the predictor domain as \code{\linkS4class{RasterBrick}} object (only exported if \code{predictor = "all"} or \code{predictor = "r_predictor"})\cr
#' \tab \code{rsq_predictor} \tab RasterBrick. Coefficient of determination between the base point and each pixel of the predictor domain as \code{\linkS4class{RasterBrick}} object (only exported if \code{predictor = "all"} or \code{predictor = "rsq_predictor"})\cr
#' \tab \code{rsq_sums_predictor} \tab RasterBrick. Sums of correlation coefficients between the base point and each pixel of the predictor domain as \code{\linkS4class{RasterBrick}} object (only exported if \code{predictor = "all"} or \code{predictor = "rsq_sums_predictor"})\cr
#' \tab \code{int_predictor} \tab RasterBrick. Intercept of the regression equation for each pixel of the predictor domain as \code{\linkS4class{RasterBrick}} object (only exported if \code{predictor = "all"} or \code{predictor = "int_predictor"})\cr
#' \tab \code{slp_predictor} \tab RasterBrick. Slope of the regression equation for each pixel of the predictor domain as \code{\linkS4class{RasterBrick}} object (only exported if \code{predictor = "all"} or \code{predictor = "slp_predictor"})\cr
#' \tab \code{r_predictor} \tab RasterBrick. The significance (p-value) of the the regression equation for each pixel of the predictor domain as \code{\linkS4class{RasterBrick}} object (only exported if \code{predictor = "all"} or \code{predictor = "p_predictor"})\cr
#' }
#' 
#' @details 
#' 
#' @author Federico Filipponi
#' 
#' @references
#' 
#' van den Dool H.M., Saha S., Johansson A. (2000). Empirical Orthogonal Teleconnections. 
#' Journal of Climate, Volume 13, Issue 8, pp. 1421-1435. 
#' \url{http://journals.ametsoc.org/doi/abs/10.1175/1520-0442%282000%29013%3C1421%3AEOT%3E2.0.CO%3B2
#' }
#' 
#' @keywords EOT time series analysis
#' 
#' @seealso \code{\link[remote]{eot}}, \code{\link[rtsa]{rtsa.eof}}, \code{\link[rtsa]{rtsa.gapfill}}, \code{\linkS4class{EOTstack}}
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
#' rmk <- pacificSST[[1]] # create raster mask
#' values(rmk) <- 1 # set raster mask values
#' rmk[which(is.na(getValues(pacificSST[[1]])))] <- 0 # set raster mask values
#'
#' ## compute EOT
#' # compute the first 2 EOTs
#' eot_result <- rtsa.eot(rasterts=rasterts, rastermask=rmk, nu=2)
#' # compute the first 2 EOTs and export only the compontents 'r_predict' and 'p_predictor'
#' eot_rp <- rtsa.eot(rasterts=rasterts, rastermask="compute", nu=2, predictor=c("r_predictor", "p_predictor"))
#' # compute the first 2 EOTs using the index of agreement
#' eot_ioa <- rtsa.eot(rasterts=rasterts, rastermask=rmk, nu=2, type="ioa", predictor="all", verbose=T)
#' }
#' 
#' @import raster
#' @import rts
#' @importFrom remote eot
#' @importFrom remote EotCycle
#' @importFrom xts as.xts
#' @importFrom methods new
#' @importFrom stats time
#' @export

# define rsta.eot function
rtsa.eot <- function(rasterts, rastermask=NULL, nu=NULL, gapfill="none", predictor="all", standardised=TRUE, reduce.both=FALSE, type="rsq", verbose=FALSE){
  
  # require the 'remote' package
  if(!requireNamespace("remote", quietly = TRUE))
    stop("Package 'remote' is required to run this function.\nYou can install using the command:\ninstall.packages('remote', dependencies=TRUE)")

  # check if input file is an object of class 'RasterStackTS', 'RasterBrickTS'
  if(!(class(rasterts) %in% c("RasterStackTS", "RasterBrickTS")))
    stop("'rasterts' argument must be an object of class 'RasterStackTS' or 'RasterBrickTS'.\nUse 'rts()' function to generate 'rasterts' input")
  
  # set number of EOTs to return
  if(is.null(nu)){
    nu <- as.integer(nlayers(rasterts))
    warning(paste("Number of EOTs to return not set using the 'nu' argument.\nIt is set by default to number of input layers: ", as.integer(rasterts@raster@data@nlayers), sep=""))
  } else {
    if(!(class(nu) %in% c("integer", "numeric"))){
      stop("'nu' argument must be numeric")
    } else {
      if(nu > as.integer(nlayers(rasterts)) | nu < 2){
        nu <- as.integer(nlayers(rasterts))
        warning(paste("Number of EOTs to return is higher than the number of input layers.\nIt is set by default to number of input layers: ", as.integer(rasterts@raster@data@nlayers), sep=""))
      }
    }
  }
  if(verbose){
    message("Number of EOTs to compute is: ", nu)
  }
  
  # check 'predictor' argument
  if(!as.logical(min(predictor %in% c("all", "r_predictor", "rsq_predictor", "rsq_sums_predictor", "int_predictor", "slp_predictor", "p_predictor")))){
    stop("'predictor' argument should contain one or more of the following: 'all', 'r_predictor', 'rsq_predictor', 'rsq_sums_predictor', 'int_predictor', 'slp_predictor', 'p_predictor'")
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
            warning("'rastermask' does not contain valid pixels for masking purpose (mask pixel values equal to 1).\nNo mask will be used to clip the input raster time series.\nSet 'rastermask' argument to 'compute' to generate it from input rasterts")
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

  # check if there are NAs in the masked dataset
  na_check <- anyNA(matrice)
  if(na_check){
    if(gapfill == "none"){
      warning("Raster time series still contain NA values after masking.\nChange input raster mask or consider the use of the available options\nfor 'gapfill' argument to select a gap-filling method before the EOT computation.\nEOT computation does not support gappy raster time series")
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
    } else {
      stop("Input rasterts contains missing values and a mask raster is not supplied.\nSet 'rastermask' argument to 'compute' to generate it from input rasterts")
    }
  } else {
    if(verbose){
      warning("Raster time series does not contain NA values after masking: OK.\nGap-filling will be not performed before the EOT computation") 
    }
  }

  # compute EOT using 'remote' package
  if(verbose){
    message(paste("EOT analysis using method '", type, "' started at: ", Sys.time(), sep=""))
    message("It may take a long time ...")
  }
  ptm <- proc.time()
  eotresult <- remote::eot(rasterts@raster, n=nu, standardised=standardised, reduce.both=reduce.both, type=type, verbose=verbose)
  if(verbose){
    message(paste(c("EOT computation ended at: "), Sys.time(), sep=""))
    message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))
  }
  
  # extract results
  # extract EOT explained variance by each mode
  eot_explained_variance <- as.vector(rep(as.double(NA), nu))
  eot_explained_variance[1] <- as.double((eotresult@modes[[1]]@cum_exp_var)*100.0)
  for(n in 2:nu){
    eot_explained_variance[n] <- as.double((eotresult@modes[[n]]@cum_exp_var - eotresult@modes[[n-1]]@cum_exp_var)*100.0)
  }
  
  # extract EOT total variance
  #eot_total_variance <- as.double(eotresult@modes[[nu]]@cum_exp_var * 100.0)
  
  # EOT Expansion Coefficient
  # create empty matrix
  eot.ec <- matrix(data=as.double(NA), nrow=length(eotresult@modes[[1]]@eot), ncol=nu)
  # fill in matrix with results
  for(l in 1:nu){
    eot.ec[,l] <- as.vector(eotresult@modes[[l]]@eot)
  }
  # change column names
  cc <- as.vector(1:nu)
  colnames(eot.ec) <- paste(c("EOT_"), sprintf("%003d", cc), sep="")
  # convert to xts object
  eot.ec <- as.xts(eot.ec, time(rasterts@time))
  
  # EOT coordinates of the base point
  eot.cbp <- matrix(data=as.double(NA), nrow=nu, ncol=2)
  colnames(eot.cbp) <- c("x", "y")
  row.names(eot.cbp) <- paste(c("EOT_"), sprintf("%003d", cc), sep="")
  for(l in 1:nu){
    eot.cbp[l,] <- as.vector(eotresult@modes[[l]]@coords_bp)
  }
  
  # create output object
  # assemble results in a object of class 'EOTstack'
  eotreturn <- new("EOTstack")
  eotreturn@eot <- eot.ec
  eotreturn@total_variance <- as.double(eotresult@modes[[nu]]@cum_exp_var * 100.0)
  eotreturn@explained_variance <- eot_explained_variance
  eotreturn@coords_bp <- eot.cbp
  
  # create output for raster object if selected with the argument 'predictor'
  
  # create empty raster stack dataset
  eot_dataset_void <- rasterts@raster[[1:nu]]
  values(eot_dataset_void) <- as.vector(rep(NA, ncell(rasterts)*nu))
  
  # mask eot r_predictor result
  if(as.logical(sum(predictor %in% c("all", "r_predictor")))){
    eot_r_predictor <- eot_dataset_void
    names(eot_r_predictor) <- paste(c("EOT_r_predictor_"), sprintf("%003d", 1:nu), sep="")
    for(l in 1:nu){
      eot_r_predictor[[l]][na_index_mask] <- eotresult@modes[[l]]@r_predictor[na_index_mask]
    }
    eotreturn@r_predictor <- eot_r_predictor
  }
  
  # mask eot rsq_predictor result
  if(as.logical(sum(predictor %in% c("all", "rsq_predictor")))){
    eot_rsq_predictor <- eot_dataset_void
    names(eot_rsq_predictor) <- paste(c("EOT_rsq_predictor_"), sprintf("%003d", 1:nu), sep="")
    for(l in 1:nu){
      eot_rsq_predictor[[l]][na_index_mask] <- eotresult@modes[[l]]@rsq_predictor[na_index_mask]
    }
    eotreturn@rsq_predictor <- eot_rsq_predictor
  }
  
  # mask eot rsq_sums_predictor result
  if(as.logical(sum(predictor %in% c("all", "rsq_sums_predictor")))){
    eot_rsq_sums_predictor <- eot_dataset_void
    names(eot_rsq_sums_predictor) <- paste(c("EOT_rsq_sums_predictor_"), sprintf("%003d", 1:nu), sep="")
    for(l in 1:nu){
      eot_rsq_sums_predictor[[l]][na_index_mask] <- eotresult@modes[[l]]@rsq_sums_predictor[na_index_mask]
    }
    eotreturn@rsq_sums_predictor <- eot_rsq_sums_predictor
  }
  
  # mask eot int_predictor result
  if(as.logical(sum(predictor %in% c("all", "int_predictor")))){
    eot_int_predictor <- eot_dataset_void
    names(eot_int_predictor) <- paste(c("EOT_int_predictor_"), sprintf("%003d", 1:nu), sep="")
    for(l in 1:nu){
      eot_int_predictor[[l]][na_index_mask] <- eotresult@modes[[l]]@int_predictor[na_index_mask]
    }
    eotreturn@int_predictor <- eot_int_predictor
  }

  # mask eot slp_predictor result
  if(as.logical(sum(predictor %in% c("all", "slp_predictor")))){
    eot_slp_predictor <- eot_dataset_void
    names(eot_slp_predictor) <- paste(c("EOT_slp_predictor_"), sprintf("%003d", 1:nu), sep="")
    for(l in 1:nu){
      eot_slp_predictor[[l]][na_index_mask] <- eotresult@modes[[l]]@slp_predictor[na_index_mask]
    }
    eotreturn@slp_predictor <- eot_slp_predictor
  }

  # mask eot p_predictor result
  if(as.logical(sum(predictor %in% c("all", "p_predictor")))){
    eot_p_predictor <- eot_dataset_void
    names(eot_p_predictor) <- paste(c("EOT_p_predictor_"), sprintf("%003d", 1:nu), sep="")
    for(l in 1:nu){
      eot_p_predictor[[l]][na_index_mask] <- eotresult@modes[[l]]@p_predictor[na_index_mask]
    }
    eotreturn@p_predictor <- eot_p_predictor
  }
  
  # mask eot resid_predictor result ### need one raster stack for each EOT. Too memory consuming!
  # eot_resid_predictor <- rasterts@raster
  # names(eot_resid_predictor) <- paste(c("EOT_resid_predictor_"), sprintf("%003d", 1:nlayers(rasterts)), sep="")
  # for(l in 1:nlayers(rasterts)){
  #   eot_resid_predictor[na_index_mask] <- eotresult@modes[[l]]@resid_predictor[na_index_mask]
  # }

  # return function result
  return(eotreturn)
  
  # stop cluster
  on.exit(stopCluster(cl))
}
