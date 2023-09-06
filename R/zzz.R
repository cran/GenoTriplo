.onAttach <- function(libname, pkgname) {
  packageStartupMessage("+-+-+-+-+-+-+-+-+-+-+")
  packageStartupMessage("|G|e|n|o|T|r|i|p|l|o|")
  packageStartupMessage("+-+-+-+-+-+-+-+-+-+-+")
  package_citation <- "Roche et al. (in prep)"
  packageStartupMessage("Thank you for using GenoTriplo!")
  packageStartupMessage("To acknowledge our work, please cite the package:")
  packageStartupMessage(package_citation)
}
