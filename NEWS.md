rtsa v0.3 (Release date: 2019-10-16)
==============

Changes:

* Added mask slot in 'EOFstack' object class (thanks RobelTakele for the fix)
* Added mask slot in 'EOTstack' object class
* Added mask raster name in in 'STDstack' object class
* Bug fix in function 'rtsa.eof' when computing only masked pixel values

rtsa v0.2 (Release date: 2018-01-25)
==============

Changes:

* Added function 'rtsa.eot' to compute Empirical Orthogonal Teleconnections
* Added function 'rtsa.stl' to compute Seasonal Trend Decomposition using Loess
* Added function 'rtsa.seas' to compute X-11 and X-13-ARIMA seasonal adjustment
* Added function 'rtsa.mk' to compute Mann-Kendall trend test

rtsa v0.1 (Release date: 2017-09-04)
==============

Changes:

* Added functions 'rtsa.eof' and 'rtsa.scaleEOF' to compute Empirical Orthogonal Functions
* Added function 'rtsa.gapfill' to perform gapfilling of missing data
