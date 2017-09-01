# title         : Methods for the 'rtsa' package
# Date          : Sep 2017
# Version       : 0.1
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>

#######################################################################

# create function to read number of layers from RasterBrickTS object
setMethod('nlayers', signature(x='RasterBrickTS'),
          function(x){
            return(x@raster@data@nlayers)
  }
)

# create function to read number of columns from RasterBrickTS object
setMethod('ncol', signature(x='RasterBrickTS'),
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

# create function to read number of cells from RasterBrickTS object
setMethod('ncell', signature(x='RasterBrickTS'),
          function(x){
            return(as.integer(as.integer(x@raster@ncols) * as.integer(x@raster@nrows)))
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
            return(x@raster@data@values)
            }
)

# # alternative method for ncell function
# # create function to read number of cells from RasterBrickTS object
# setMethod('ncell', signature(x='RasterBrickTS'),
#           function(x){
#             return(nrow(x) * ncol(x))
#           }
# )

### note

#nlayers()
#if raster is inMemory the number of layer can be retrieved using the method above
#if raster is not in memory number of layer can be retrieved using:
#length(rasterts@raster@layers)
# as an alternative the number of layers can be derived from the length of dates in the rts object
#length(rasterts@time[,1])

#ncell()
#by default rts package return value '1' when using ncell()