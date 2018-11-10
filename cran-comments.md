Hello, 

I am resubmitting the package after fixing the 1 NOTE about the non-standard directory docs and file cran-comments.md present in the built package (added to .Rbuildignore). Please, ignore the note about misspelled words - they are all correct.

# Output from check_rhub:

The downloaded binary packages are in
	/var/folders/nb/_b5mkf753p3c86s9nrpk33gc002289/T//RtmpszhHI7/downloaded_packages
✔  checking for file ‘/Users/vmitov/Documents/Bio/Projects/PCMBaseR/DESCRIPTION’ ...
─  preparing ‘PCMBase’: (681ms)
✔  checking DESCRIPTION meta-information ...
─  installing the package to build vignettes (375ms)
✔  creating vignettes (8.2s)
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘PCMBase_1.2.7.tar.gz’
   Warning: invalid uid value replaced by that for user 'nobody'
   
─  Uploading package
─  Preparing build, see status at
   http://builder.r-hub.io/status/PCMBase_1.2.7.tar.gz-30f36f36bb2640ac899dc960a21c111d
   http://builder.r-hub.io/status/PCMBase_1.2.7.tar.gz-f2ccbd3838c2402bb984becf1e1bd6f0
   http://builder.r-hub.io/status/PCMBase_1.2.7.tar.gz-405a62181dd847fc9e3dce58b70b6ab9
─  Build started
─  Creating new user
─  Downloading and unpacking package file
─  Querying package dependencies
─  Installing package dependencies
─  Running R CMD check
   setting _R_CHECK_FORCE_SUGGESTS_ to false
   setting R_COMPILE_AND_INSTALL_PACKAGES to never
   setting R_REMOTES_STANDALONE to true
   setting R_REMOTES_NO_ERRORS_FROM_WARNINGS to true
   setting _R_CHECK_FORCE_SUGGESTS_ to true
   setting _R_CHECK_CRAN_INCOMING_USE_ASPELL_ to true
─  using log directory 'C:/Users/USERSxJciIZGHR/PCMBase.Rcheck'
─  using R Under development (unstable) (2018-09-27 r75377)
─  using platform: x86_64-w64-mingw32 (64-bit)
─  using session charset: ISO8859-1 (1.4s)
─  using option '--as-cran'
✔  checking for file 'PCMBase/DESCRIPTION'
─  checking extension type ... Package
─  this is package 'PCMBase' version '1.2.7' (801ms)
─  package encoding: UTF-8
N  checking CRAN incoming feasibility
   Maintainer: 'Venelin Mitov <vmitov@gmail.com>'
   
   New submission
   
   Possibly mis-spelled words in DESCRIPTION:
     OU (18:56, 26:62, 26:69, 27:13)
     Ornstein (18:36)
     PCMBase (21:17, 28:75)
     PCMFit (30:28)
     Uhlenbeck (18:45)
✔  checking package namespace information
✔  checking package dependencies
✔  checking if this is a source package
✔  checking if there is a namespace (1.6s)
✔  checking for executable files
✔  checking for hidden files and directories
✔  checking for portable file names
✔  checking serialization versions (808ms)
✔  checking whether package 'PCMBase' can be installed
✔  checking installed package size
✔  checking package directory
✔  checking for future file timestanps (805ms)
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information
✔  checking top-level files
✔  checking for left-over files (803ms)
✔  checking index information
✔  checking package subdirectories
✔  checking R files for non-ASCII characters
✔  checking R files for syntax errors (1.6s)
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly
✔  checking whether the namespace can be loaded with stated dependencies (807ms)
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path
✔  checking use of S3 registration
✔  checking dependencies in R code (800ms)
✔  checking S3 generic/method consistency
✔  checking replacement functions
✔  checking foreign function calls
✔  checking R code for possible problems (11.2s)
✔  checking Rd files (803ms)
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references (799ms)
✔  checking for missing documentation entries (2.4s)
✔  checking for code/documentation mismatches (8.8s)
✔  checking Rd \usage sections (4.8s)
✔  checking Rd contents (801ms)
✔  checking for unstated dependencies in examples
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
─  checking examples ... NONE
✔  checking for unstated dependencies in vignettes (1.6s)
✔  checking package vignettes in 'inst/doc'
✔  checking re-building of vignette outputs (6.4s)
✔  checking PDF version of manual (20.1s)
   
─  Done with R CMD check
─  Cleaning up files and user
    

── PCMBase 1.2.7: NOTE

  Build ID:   PCMBase_1.2.7.tar.gz-30f36f36bb2640ac899dc960a21c111d
  Platform:   Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  Submitted:  4.4s ago
  Build time: 3m 57.2s

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Venelin Mitov <vmitov@gmail.com>'
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    OU (18:56, 26:62, 26:69, 27:13)
    Ornstein (18:36)
    PCMBase (21:17, 28:75)
    PCMFit (30:28)
    Uhlenbeck (18:45)

0 errors ✔ | 0 warnings ✔ | 1 note ✖

── PCMBase 1.2.7: IN-PROGRESS

  Build ID:   PCMBase_1.2.7.tar.gz-f2ccbd3838c2402bb984becf1e1bd6f0
  Platform:   Ubuntu Linux 16.04 LTS, R-release, GCC
  Submitted:  4.4s ago


── PCMBase 1.2.7: IN-PROGRESS

  Build ID:   PCMBase_1.2.7.tar.gz-405a62181dd847fc9e3dce58b70b6ab9
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  4.4s ago

Venelin Mitov: All of the words are correclty spelled.

── PCMBase 1.2.7: NOTE

  Build ID:   PCMBase_1.2.7.tar.gz-3c1a4aa5287e4677b5b12df498e01c05
  Platform:   Ubuntu Linux 16.04 LTS, R-release, GCC
  Submitted:  6.5s ago
  
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Venelin Mitov <vmitov@gmail.com>’

New submission

Possibly mis-spelled words in DESCRIPTION:
  Ornstein (18:36)
  OU (18:56, 26:62, 26:69, 27:13)
  PCMBase (21:17, 28:75)
  PCMFit (30:28)
  Uhlenbeck (18:45)
  

