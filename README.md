# rtsa
## R package for Raster Time Series Analysis

The package provides a collection of analytics to perform spatio-temporal analysis from raster time series. 
It acts as a front-end to already available functions in various R packages, specifically designed to handle geographic datasets provided as raster time series. The available functions within the package allow the direct input of raster time series to extract concise and comprehensive information.
Since some techniques for spatio-temporal analysis can not deal with missing values raster time series, a selection of gap-filling methods are provided.

### Main features:

* use of raster time series with explicit temporal dimension
* use of raster time series as direct input to functions
* use of a raster mask to select the region of interest and reduce memory loads
* parallel processing using multiple CPUs

### Supported gap-filling methods:

* DINEOF
* linear interpolation
* spline interpolation
* stine interpolation

### Currently, the following analytical methods are supported:

* Empirical Orthogonal Function
* Empirical Orthogonal Teleconnections
* Seasonal Trend Decomposition using Loess
* X-11
* X-13-ARIMA seasonal adjustment
* Mann-Kendall trend test

### The following analytical methods will be supported soon:

* Self Organizing Maps

### Installation

**To load** (using `devtools`):
```r
library(devtools)
install_github("marchtaylor/sinkr")
install_github("ffilipponi/rtsa")
```

### Authors

* Filipponi Federico

### License

Licensed under the GNU General Public License, Version 3.0: https://www.gnu.org/licenses/gpl-3.0.html

### Funding

The ECOPOTENTIAL project has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 641762
