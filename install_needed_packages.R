# Code in this directory will help to install required packages from CRAN or Bioconductor repositories


#' For Cran packages
#pkgs is a vector of package names that are required from CRAN
install_cran_pkgs <- function(pkgs){
  # installed.packages()[,"Package"] will list all the packages that are installed
  # The whole code will create vector not_installed with missing packages.
  missing_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  # if statement checks if missing_pkgs have any elements in it.
  # if there are it proceeds inside the function for installation.
  if(length(missing_pkgs)){
    install.packages(missing_pkgs, dependencies = TRUE) 
    # you may get a warning for ewas here, just ignore that
    # as ewastools is available from github
    if("ewastools" %in% pkgs & "ewastools" %in% missing_pkgs){
      # Package source : https://github.com/hhhh5/ewastools
      devtools::install_github("hhhh5/ewastools@master")
    }
  }else{
    # missing_pkgs vector is empty
    message("Cran packages you need are installed")
  }
}

#' For bioconductor packages
# The code executes exactly as above except it uses install function of
# BiocManager to install the packages.
install_bioconductor_pkgs <- function(pkgs){
  bioc_missing_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(bioc_missing_pkgs)){
    BiocManager::install(c(bioc_missing_pkgs))
  }else{
    message("Bioconductor packages you need are installed")
  }
}


#' Run test to check if all packages are installed
#' If not, figure out the problem and make sure all packages are istalled
check_installed <- function(pkgs){
  not_installed <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(not_installed)){
    message("Package(s) missing from list of installed packages: ", not_installed)
    message("Please install missing packages manually and sort any issues")
    message("Hopefully there are not many .... GOOD LUCK !!!")
  }else{
    message("All Required Packages are installed")
    message("Awesome lets have fun !!!!!")
  }
}




