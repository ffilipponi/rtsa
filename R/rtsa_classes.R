# title         : Classes for the 'rtsa' package
# Date          : Sep 2017
# Version       : 0.1
# Licence       : GPL v3
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>

#######################################################################

setClass("EOFstack",
         representation(eof="RasterBrick",
                        expansion_coefficients="xts",
                        total_variance="numeric",
                        explained_variance="numeric",
                        center="RasterLayer",
                        scale="RasterLayer"),
         prototype=prototype(eof=NULL, expansion_coefficients=NULL, total_variance=NULL, explained_variance=NULL, center=NULL, scale=NULL)
         )
