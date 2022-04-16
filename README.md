# PlexedPiper

<!-- badges: start -->
[![R-CMD-check](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/workflows/R-CMD-check/badge.svg)](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/actions)
<!-- badges: end -->


R package used at PNNL for processing isobaric labeling (e.g. TMT) proteomics
data. The key inputs are:

* MS/MS identifications from the [MS-GF+ search engine](https://github.com/MSGFPlus/msgfplus)
* Reporter ion intensities extracted using [MASIC](https://github.com/pnnl-comp-mass-spec/MASIC)
* Tables outlining study design
   + table linking dataset to plexes
   + table linking reporter channels with sample names
   + table identifying reference within each plex

## Updates

Updates to this package are detailed in the [NEWS.md](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/blob/master/NEWS.md) file and in the [releases](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/releases).

## Website

The [PlexedPiper pkgdown website](https://pnnl-comp-mass-spec.github.io/PlexedPiper/) contains pre-built vignettes as well as updates and documentation.

## R Installation and Usage

```R
if(!require("remotes", quietly = T)) install.packages("remotes")
remotes::install_github("PNNL-Comp-Mass-Spec/PlexedPiper", build_vignettes = TRUE)
library(PlexedPiper)
vignette("tmt_pipeline_v1")
```

### Example Data

A companion R package with test data based on the MoTrPAC pilot study is available
here

* [PlexedPiperTestData](https://github.com/vladpetyuk/PlexedPiperTestData)

## Docker/Linux installation

PlexedPiper can be run within a [Docker Container](https://www.docker.com/resources/what-container)

* This example `Dockerfile` shows the required system libraries, starting with the base [rocker/TidyVerse](https://hub.docker.com/r/rocker/tidyverse/dockerfile) image

```Dockerfile
FROM rocker/tidyverse:3.6.1
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 unixodbc \
 unixodbc-dev \
 freetds-dev \
 freetds-bin \
 tdsodbc \
 libcurl4-openssl-dev \
 libxml2-dev \
 libnetcdf-dev \
 libssl-dev
RUN R -e 'remotes::install_github("PNNL-Comp-Mass-Spec/PlexedPiper", build_vignettes = TRUE)'
```

## MacOS installation

On MacOS, install [Homebrew](https://brew.sh/), then use

```Shell
brew install unixodbc
brew install freetds
```
Note, the `--with-unixodbc` option in freetds installation is deprecated.

Create `~/.odbcinst.ini` file and add
```INI
[FreeTDS]
Driver = /usr/local/lib/libtdsodbc.so
```
If your location of `libtdsodbc.so` differs, use the proper location.

### Installation Tips

If within PNNL network there may be an error associated with `mount_smbfs`. This happens due to network access credentials. Options are either to wait or proactively access one of the PNNL servers. For example try mounting one of the public directories from the terminal window. Enter your network password once requested. 
`mount -t smbfs //protoapps/DataPkgs/Public/ ~/temp_msms_results`
Then compilation of the vignettes that imply access to PNNL DMS should proceed smoothly.


Original Location
====
The original location is on the [vladpetyuk](https://github.com/vladpetyuk) account, repo [PlexedPiper](https://github.com/vladpetyuk/PlexedPiper).




