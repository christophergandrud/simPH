## R CMD check results

0 errors | 1 warning | 1 note

* checking sizes of PDF files under 'inst/doc' ... WARNING
  'gs+qpdf' made some significant size reductions:
     compacted 'simPH-overview.pdf' from 677Kb to 305Kb
  consider running tools::compactPDF(gs_quality = "ebook") on these files,
  or build the source package with --compact-vignettes=both
  
* checking for future file timestamps ... NOTE
  unable to verify current time
  
### Comments on check results:
* The WARNING about PDF sizes will be addressed before final submission by building the package with --compact-vignettes=both
* The NOTE about future file timestamps is due to system time verification issues and does not affect package functionality
* The NOTE about "Non-standard file/directory found at top level: 'cran-comments.md'" is expected and can be ignored as this file is for CRAN submission purposes

## Release summary

This release (1.3.14):

* Improves documentation to meet new CRAN requirements
* Updates class() comparisons to use inherits() instead
* Fixes issues with LaTeX equations in documentation
* Addresses various documentation formatting issues

## Test environments

* local macOS install, R 4.3.1
* win-builder (devel and release) - will be tested before CRAN submission
* R-hub platforms - will be tested before CRAN submission

## Downstream dependencies

There are no downstream dependencies for this package.
