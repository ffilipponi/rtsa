# title         : Classes for the 'rtsa' package
# Date          : Sep 2017
# Version       : 0.1
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
         representation(eof="RasterBrick",
                        expansion_coefficients="xts",
                        total_variance="numeric",
                        explained_variance="numeric",
                        center="RasterLayer",
                        scale="RasterLayer"),
         prototype=prototype(eof=NULL, expansion_coefficients=NULL, total_variance=NULL, explained_variance=NULL, center=NULL, scale=NULL)
         )
