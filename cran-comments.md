## R CMD check results

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking top-level files ... NOTE
  Non-standard file/directory found at top level:
    ‘cran-comments.md’

0 errors ✔ | 0 warnings ✔ | 2 notes ✖
  
  Non-standard file/directory found at top level:
  'cran-comments.md'
  
### Comments on check results:
* The NOTE about future file timestamps is due to system time verification issues and does not affect package functionality
* The NOTE about "Non-standard file/directory found at top level: 'cran-comments.md'" is expected and can be ignored as this file is for CRAN submission purposes

## Release summary

This release (1.3.14):

* Improves documentation to meet new CRAN requirements
* Updates class() comparisons to use inherits() instead
* Fixes issues with LaTeX equations in documentation
* Updates the CITATION file to use modern syntax
* Fixes redirected URLs in documentation
* Updates package anchor links in Rd files
* Adds vignette compression
* Addresses various documentation formatting issues

## Test environments

* local macOS install, R 4.3.1
* win-builder (devel and release) - will be tested before CRAN submission
* R-hub platforms - will be tested before CRAN submission

## Downstream dependencies

There are no downstream dependencies for this package.
