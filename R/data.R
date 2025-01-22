#' Converter for TMT Reporter Channels
#'
#' Converter from human-readable notation of TMT reporter channels to
#' m/z-type (reported by MASIC).
#'
#' @format A data frame with 10 rows and 2 variables:
#' \describe{
#'   \item{ReporterName}{human-readable notation}
#'   \item{ReporterIon}{reporter channel names extracted by MASIC}
#'   ...
#' }
#' @source \url{https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-mass-spectrometry-analysis/protein-quantitation-mass-spectrometry/tandem-mass-tag-systems.html}
"reporter_converter"


#' Converter for new TMT masses
#'
#' If any reporter ions have new reported masses, this list will convert them to the old
#'    masses.
#'
#' @format A list with names containing new reported masses and values containing old 
#'    reported massess
#' @source \url{https://assets.thermofisher.com/TFS-Assets%2FBID%2Fposters%2Fexpanding-tmtpro-reagents-32-plex-orbitrap-platforms-poster.pdf}
"reporter_converter_new_masses"
