# title         : Classes for the 'rtsa' package
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>

#######################################################################
#
#' Class "EOFstack"
#'  
#' @name EOFstack
#' @description Result from \code{\link[rtsa]{rtsa.eof}} stored as a \code{\linkS4class{EOFstack}} object.
#' Slots for \code{EOFstack} objects include: (1) the spatial representation of eof modes; 
#' (2) the eof expansion coefficients; (3) the dataset total variance;
#' (4) the explained variance from each eof; (5) the raster of center values; 
#' (6) the raster of scale values.
#' @docType class
#' @section Objects from the class: Objects are created by calls to \code{\link[rtsa]{rtsa.eof}}.
#' @seealso \code{\link[rtsa]{rtsa.eof}}

#' @exportClass EOFstack

setClass("EOFstack",
         representation(eof.modes="RasterBrick",
                        expansion_coefficients="xts",
                        total_variance="numeric",
                        explained_variance="numeric",
                        center="RasterLayer",
                        scale="RasterLayer"),
         prototype=prototype(eof.modes=NULL, expansion_coefficients=NULL, total_variance=NULL, explained_variance=NULL, center=NULL, scale=NULL)
         )

setClass("EOTstack",
         representation(eot="xts",
                        total_variance="numeric",
                        explained_variance="numeric",
                        coords_bp="matrix",
                        r_predictor="RasterBrick",
                        rsq_predictor="RasterBrick",
                        rsq_sums_predictor="RasterBrick",
                        int_predictor="RasterBrick",
                        slp_predictor="RasterBrick",
                        p_predictor="RasterBrick"),
         prototype=prototype(eot=NULL, total_variance=NULL, explained_variance=NULL, coords_bp=NULL, r_predictor=NULL, rsq_predictor=NULL, rsq_sums_predictor=NULL, int_predictor=NULL, slp_predictor=NULL, p_predictor=NULL)
)

setClass("STDstack",
         representation(std="character",
                        mask="RasterLayer",
                        seasonal_amplitude="RasterLayer",
                        seasonal_amplitude_stdev="RasterLayer",
                        trend_slope="RasterLayer",
                        remainder_stdev="RasterLayer",
                        rts="RasterBrickTS",
                        seasonality="RasterBrickTS",
                        trend="RasterBrickTS",
                        seasonaladjtrend="RasterBrickTS",
                        remainder="RasterBrickTS"),
         prototype=prototype(std=NULL, rts=NULL, seasonality=NULL, trend=NULL, seasonaladjtrend=NULL, remainder=NULL)
)

setClass("MKstack",
         representation(tau="RasterLayer",
                        pvalue="RasterLayer",
                        score="RasterLayer",
                        variance="RasterLayer"),
         prototype=prototype(tau=NULL, pvalue=NULL, score=NULL, variance=NULL)
)

setClass("SOMstack",
         representation(som="RasterBrick",
                        bmu="xts"),
         prototype=prototype(som=NULL, bmu=NULL)
)
