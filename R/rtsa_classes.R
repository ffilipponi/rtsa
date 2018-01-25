# title         : Classes for the 'rtsa' package
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
#
#######################################################################
#
#' Class "EOFstack"
#'  
#' @name EOFstack-class
#' @description Result from \code{\link[rtsa]{rtsa.eof}} stored as a \code{\linkS4class{EOFstack}} object.
#' Slots for \code{EOFstack} objects include:
#' \itemize{
#' \item \emph{eof.modes}: the spatial representation of eof modes;
#' \item \emph{expansion_coefficients}: the eof expansion coefficients;
#' \item \emph{total_variance}: the dataset total variance;
#' \item \emph{explained_variance}: the explained variance by each eof mode;
#' \item \emph{center}: the raster of center values;
#' \item \emph{scale}: the raster of scale values.
#' }
#' @docType class
#' @section Objects from the class: Objects are created by calls to \code{\link[rtsa]{rtsa.eof}}.
#' @seealso \code{\link[rtsa]{rtsa.eof}}
#' 
#' @exportClass EOFstack
#' @rdname EOFstack-class

setClass("EOFstack",
         representation(eof.modes="RasterBrick",
                        expansion_coefficients="xts",
                        total_variance="numeric",
                        explained_variance="numeric",
                        center="RasterLayer",
                        scale="RasterLayer"),
         prototype=prototype(eof.modes=NULL, expansion_coefficients=NULL, total_variance=NULL, explained_variance=NULL, center=NULL, scale=NULL)
)

NULL

#' Class "EOTstack"
#' 
#' @name EOTstack-class
#' @description Result from \code{\link[rtsa]{rtsa.eot}} stored as a \code{\linkS4class{EOTstack}} object.
#' Slots for \code{EOTstack} objects include: 
#' \itemize{
#' \item \emph{eot}: the temporal profiles of the identified base point pixels;
#' \item \emph{total_variance}: the dataset total variance;
#' \item \emph{explained_variance}: the explained variance by each eot mode;
#' \item \emph{coords_bp}: the coordinates of the identified base point;
#' \item \emph{r_predictor}: the raster of correlation coefficients between the base point and each pixel of the predictor domain;
#' \item \emph{rsq_predictor}: the raster of coefficient of determination between the base point and each pixel of the predictor domain;
#' \item \emph{rsq_sums_predictor}: the raster of sums of coefficient of determination between the base point and each pixel of the predictor domain;
#' \item \emph{int_predictor}: the raster of intercept of the regression equation for each pixel of the predictor domain;
#' \item \emph{slp_predictor}: the raster of slope of the regression equation for each pixel of the predictor domain;
#' \item \emph{p_predictor}: the raster of significance (p-value) of the the regression equation for each pixel of the predictor domain.
#' }
#' @docType class
#' @section Objects from the class: Objects are created by calls to \code{\link[rtsa]{rtsa.eot}}.
#' @seealso \code{\link[rtsa]{rtsa.eot}}
#' 
#' @exportClass EOTstack
#' @rdname EOTstack-class
#' 
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

NULL

#' Class "STDstack"
#' 
#' @name STDstack-class
#' @description Result from \code{\link[rtsa]{rtsa.stl}} and \code{\link[rtsa]{rtsa.seas}} stored as a \code{\linkS4class{STDstack}} object.
#' Slots for \code{STDstack} objects include: 
#' \itemize{
#' \item \emph{std}: Seasonal Trend Decomposition method used;
#' \item \emph{mask}: the raster mask of valid pixels;
#' \item \emph{seasonal_amplitude}: the raster of seasonal amplitude;
#' \item \emph{seasonal_amplitude}: the raster of seasonal amplitude standard deviation;
#' \item \emph{trend_slope}: the raster of trend slope value;
#' \item \emph{remainder_stdev}: the raster of the remainder standard deviation;
#' \item \emph{rts}: input raster time series as RasterBrickTS object;
#' \item \emph{seasonality}: seasonal component as RasterBrickTS object;
#' \item \emph{trend}: trend component as RasterBrickTS object;
#' \item \emph{seasonaladjtrend}: seasonal adjusted trend component as RasterBrickTS object;
#' \item \emph{remainder}: remainder component as RasterBrickTS object.
#' }
#' @docType class
#' @section Objects from the class: Objects are created by calls to \code{\link[rtsa]{rtsa.stl}} and \code{\link[rtsa]{rtsa.seas}}.
#' @seealso \code{\link[rtsa]{rtsa.stl}}, \code{\link[rtsa]{rtsa.seas}}
#' 
#' @exportClass STDstack
#' @rdname STDstack-class
#' 
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

NULL

#' Class "MKstack"
#' 
#' @name MKstack-class
#' @description Result from \code{\link[rtsa]{rtsa.mk}} stored as a \code{\linkS4class{MKstack}} object.
#' Slots for \code{MKstack} objects include: 
#' \itemize{
#' \item \emph{tau}: the raster of Kendall tau statistic;
#' \item \emph{pvalue}: the raster of Kendall two-sided p-value;
#' \item \emph{score}: the raster of Kendall Score;
#' \item \emph{variance}: the raster of Variance of Kendall Score.
#' }
#' @docType class
#' @section Objects from the class: Objects are created by calls to \code{\link[rtsa]{rtsa.mk}}.
#' @seealso \code{\link[rtsa]{rtsa.stl}}, \code{\link[rtsa]{rtsa.seas}}, \code{\link[rtsa]{rtsa.mk}}
#' 
#' @exportClass MKstack
#' @rdname  MKstack-class
#' 
setClass("MKstack",
         representation(tau="RasterLayer",
                        pvalue="RasterLayer",
                        score="RasterLayer",
                        variance="RasterLayer"),
         prototype=prototype(tau=NULL, pvalue=NULL, score=NULL, variance=NULL)
)

NULL
 
# #' Class "SOMstack"
# #' 
# #' @name SOMstack-class
# #' @description Result from \code{\link[rtsa]{rtsa.som}} stored as a \code{\linkS4class{SOMstack}} object.
# #' Slots for \code{SOMstack} objects include: 
# #' \itemize{
# #' \item \emph{som}: the spatial representation of resulting SOMs;
# #' \item \emph{bmu}: Best Matching Unit for each input temporal observation.
# #' }
# #' @docType class
# #' @section Objects from the class: Objects are created by calls to \code{\link[rtsa]{rtsa.som}}.
# #' @seealso \code{\link[rtsa]{rtsa.som}}
# #'
# #' @exportClass SOMstack
# #' @rdname SOMstack-class
# 
# setClass("SOMstack",
#          representation(som="RasterBrick",
#                         bmu="xts"),
#          prototype=prototype(som=NULL, bmu=NULL)
# )
