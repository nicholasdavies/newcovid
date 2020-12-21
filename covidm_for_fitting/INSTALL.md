# Installation guide

Install all packages required in the script.

```R
install.packages(c("Rcpp" , "RcppGSL", "HDInterval", "ggplot2", 
                   "data.table", "socialmixr", "shiny", "lubridate", 
                   "readxl", "cowplot", "stringr", "Hmisc", "extraDistr",
                   "nloptr", "viridis", "magick"))
```

Additional C++ libraries are required.

## OS X

1. Ensure you have access to `g++` at your terminal.
2. Download and install [CMake](https://cmake.org/download/)
3. Download and install [gfortran](https://gcc.gnu.org/wiki/GFortranBinariesMacOS)
4. Download and install [nlopt](https://nlopt.readthedocs.io/en/latest/NLopt_Installation/)
   1. If `cmake . `  doesn't immediately work, you might need to add CMake to your path, using the following code in the path `PATH="/Applications/CMake.app/Contents/bin":"$PATH"` 

## Linux

1. Using your favourite package manager, install `cmake`, `gfortran`, and `libnlopt-cxx-dev`
2. If you don't have the R packages listed above installed (or if they won't install) you can install them from your favourite package manager where they're probably called `r-cran-XYZ` where `XYZ` is the name of the R package

