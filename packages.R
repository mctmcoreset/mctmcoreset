#packages.R
library("tram") 
library("mvtnorm")
library("colorspace")
library("latex2exp")
library("Matrix")
library("pracma")
library(geometry)


# List of packages
packages <- c("tidyr", "ggplot2", "egg", "bigmemory", "Matrix", "pracma",
              "tram","mvtnorm","colorspace","latex2exp","Matrix",
              "pracma","geometry")

# Check and install (if not already) and then load the package
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

