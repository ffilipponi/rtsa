# title         : Methods for the 'rtsa' package
# Date          : Jan 2018
# Version       : 0.2
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>

#######################################################################

#' @importFrom methods .hasSlot
#' @exportMethod nlayers ncol nrow ncell setValues getValues

setMethod('nlayers', signature(x='RasterBrickTS'),
          function(x){
            return(length(x@time))
  }
)

setMethod('nlayers', signature(x='RasterStackTS'),
          function(x){
            return(length(x@time))
          }
)

# create function to read number of columns from RasterBrickTS object
setMethod('ncol', signature(x='RasterBrickTS'),
          function(x){
            return(x@raster@ncols)
          }
)

# create function to read number of columns from RasterBrickTS object
setMethod('ncol', signature(x='RasterStackTS'),
          function(x){
            return(x@raster@ncols)
          }
)

# create function to read number of rows from RasterBrickTS object
setMethod('nrow', signature(x='RasterBrickTS'),
          function(x){
            return(x@raster@nrows)
          }
)

# create function to read number of columns from RasterBrickTS object
setMethod('nrow', signature(x='RasterStackTS'),
          function(x){
            return(x@raster@ncols)
          }
)

# create function to read number of cells from RasterBrickTS object
setMethod('ncell', signature(x='RasterBrickTS'),
          function(x){
            return(ncol(x) * nrow(x))
          }
)

# create function to read number of cells from RasterStackTS object
setMethod('ncell', signature(x='RasterStackTS'),
          function(x){
            return(ncol(x) * nrow(x))
          }
)

# create function to set cell values to RasterBrickTS object
setMethod('setValues', signature(x='RasterBrickTS'),
          function(x, values){
            if(length(rasterts@raster@data@values) != length(as.vector(values))){
              stop("Length of values differs from target RasterBrickTS size")
            } else {
              x@raster@data@values <- matrix(data=as.vector(values), nrow=ncell(x), ncol=nlayers(x))
              return(x)
            }
          }
)

# create function to get cell values to RasterBrickTS object
setMethod('getValues', signature(x='RasterBrickTS'),
          function(x){
            if(.hasSlot(x@raster, 'data')){
              return(x@raster@data@values)
            } else {
              stop("Input argument does not contain a 'data' slot.\nUse 'inMemory' raster object as input for the 'rts' function")
            }
            }
)

setMethod('getValues', signature(x='RasterStackTS'),
          function(x){
            if(.hasSlot(x@raster, 'data')){
              return(x@raster@data@values)
            } else {
              stop("Input argument does not contain a 'data' slot.\nUse 'inMemory' raster object as input for the 'rts' function")
            }
          }
)
