Hello, 

I am resubmitting the package after fixing the 1 NOTE that there were no examples in the Rd pages. I've now written examples to the most important functions. Additional examples are available in the package web-site and vignettes. 

Previous notes:

For running the rhub check I've used the command: 

rhub::check_for_cran(env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false", `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"))

Setting _R_CHECK_FORCE_SUGGESTS_ to false was necessary because of the package ggtree, which is listed in Suggests and is used conditionally in one function of the package. The ggtree package is hosted on Bioconductor.

There was a problem with the installation of the package on Fedora (during check_rhub()). I could find the following error in the log-file (the whole log-file is available from this link: https://builder.r-hub.io/status/original/PCMBase_1.2.7.tar.gz-405a62181dd847fc9e3dce58b70b6ab9):

Error response from daemon: No such container: PCMBase_1.2.7.tar.gz-405a62181dd847fc9e3dce58b70b6ab9-3
+

Please, ignore the note about misspelled words - they are all correct. 

# Output from check_win_devel()
* using log directory 'd:/RCompile/CRANguest/R-devel/PCMBase.Rcheck'
* using R Under development (unstable) (2018-11-08 r75566)
* using platform: x86_64-w64-mingw32 (64-bit)
* using session charset: ISO8859-1
* checking for file 'PCMBase/DESCRIPTION' ... OK
* checking extension type ... Package
* this is package 'PCMBase' version '1.2.7'
* package encoding: UTF-8
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Venelin Mitov <vmitov@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  OU (18:56, 26:62, 26:69, 27:13)
  Ornstein (18:36)
  PCMBase (21:17, 28:75)
  PCMFit (30:28)
  Uhlenbeck (18:45)
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking serialization versions ... OK
* checking whether package 'PCMBase' can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking 'build' directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* loading checks for arch 'i386'
** checking whether the package can be loaded ... OK
** checking whether the package can be loaded with stated dependencies ... OK
** checking whether the package can be unloaded cleanly ... OK
** checking whether the namespace can be loaded with stated dependencies ... OK
** checking whether the namespace can be unloaded cleanly ... OK
** checking loading without being on the library search path ... OK
** checking use of S3 registration ... OK
* loading checks for arch 'x64'
** checking whether the package can be loaded ... OK
** checking whether the package can be loaded with stated dependencies ... OK
** checking whether the package can be unloaded cleanly ... OK
** checking whether the namespace can be loaded with stated dependencies ... OK
** checking whether the namespace can be unloaded cleanly ... OK
** checking loading without being on the library search path ... OK
** checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... [19s] OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking installed files from 'inst/doc' ... OK
* checking files in 'vignettes' ... OK
* checking examples ... NONE
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in 'inst/doc' ... OK
* checking re-building of vignette outputs ... [6s] OK
* checking PDF version of manual ... OK
* DONE
Status: 1 NOTE

# Output from check_rhub:
─  Building package
─  Uploading package
─  Preparing build, see status at
   http://builder.r-hub.io/status/PCMBase_1.2.7.tar.gz-15019b7100bb41ed831205526f445267
   http://builder.r-hub.io/status/PCMBase_1.2.7.tar.gz-abde6a3fb3d04e1ab193b31d05b7755f
   http://builder.r-hub.io/status/PCMBase_1.2.7.tar.gz-7c59ac3936fa45b28eec46d345e223eb
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
   setting _R_CHECK_FORCE_SUGGESTS_ to false
   setting _R_CHECK_CRAN_INCOMING_USE_ASPELL_ to true
─  using log directory 'C:/Users/USEREXBIKcvLzr/PCMBase.Rcheck'
─  using R Under development (unstable) (2018-09-27 r75377)
─  using platform: x86_64-w64-mingw32 (64-bit)
─  using session charset: ISO8859-1 (1.6s)
─  using option '--as-cran'
✔  checking for file 'PCMBase/DESCRIPTION'
─  checking extension type ... Package
─  this is package 'PCMBase' version '1.2.7' (797ms)
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
N  checking package dependencies
   Package suggested but not available for checking: 'ggtree'
✔  checking if this is a source package (793ms)
✔  checking if there is a namespace
✔  checking for executable files
✔  checking for hidden files and directories
✔  checking for portable file names (791ms)
✔  checking serialization versions
✔  checking whether package 'PCMBase' can be installed
✔  checking installed package size
✔  checking package directory (1.6s)
✔  checking for future file timestanps
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information
✔  checking top-level files (806ms)
✔  checking for left-over files
✔  checking index information
✔  checking package subdirectories
✔  checking R files for non-ASCII characters (790ms)
✔  checking R files for syntax errors
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly (796ms)
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path
✔  checking use of S3 registration (1.6s)
✔  checking dependencies in R code
✔  checking S3 generic/method consistency
✔  checking replacement functions
✔  checking foreign function calls (789ms)
✔  checking R code for possible problems (7.1s)
✔  checking Rd files (804ms)
✔  checking Rd metadata
✔  checking Rd line widths (1.6s)
✔  checking Rd cross-references
✔  checking for missing documentation entries (2.4s)
✔  checking for code/documentation mismatches (7.9s)
✔  checking Rd \usage sections (4.7s)
✔  checking Rd contents
✔  checking for unstated dependencies in examples (2.4s)
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
─  checking examples ... NONE
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes in 'inst/doc'
✔  checking re-building of vignette outputs (5.6s)
✔  checking PDF version of manual (21.4s)
   
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

── PCMBase 1.2.7: NOTE

  Build ID:   PCMBase_1.2.7.tar.gz-f2ccbd3838c2402bb984becf1e1bd6f0
  Platform:   Ubuntu Linux 16.04 LTS, R-release, GCC
  Submitted:  4.4s ago


── PCMBase 1.2.7: ERROR 

  Build ID:   PCMBase_1.2.7.tar.gz-405a62181dd847fc9e3dce58b70b6ab9
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  4.4s ago


