.onAttach <- function(libname, pkgname) {
  packageStartupMessage("+-+-+-+-+-+-+-+-+-+-+")
  packageStartupMessage("|G|e|n|o|T|r|i|p|l|o|")
  packageStartupMessage("+-+-+-+-+-+-+-+-+-+-+")
  package_citation = "Roche et al. (2024). GenoTriplo: A SNP genotype calling method for triploids."
  doi = "https://doi.org/10.1371/journal.pcbi.1012483"
  packageStartupMessage("Thank you for using GenoTriplo!")
  packageStartupMessage("To acknowledge our work, please cite the package:")
  packageStartupMessage(package_citation)
  packageStartupMessage(doi)
}
