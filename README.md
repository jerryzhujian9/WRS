# WRS #

A package of R.R. Wilcox' robust statistics functions.
For more information, see [http://dornsife.usc.edu/labs/rwilcox/software/](http://dornsife.usc.edu/labs/rwilcox/software/).

Jerry: forked at Tue, Jan 16 2018
http://dornsife.usc.edu/assets/sites/239/docs/Rallfun-v34.txt

merged additional C functions from http://dornsife.usc.edu/assets/sites/239/docs/WRSC.txt
    lintest_C()
    scorsubMC()
    scorci_C()

## Installation ##

    install.packages( c("Rcpp","RcppArmadillo") )
    devtools::install_github('jerryzhujian9/WRScpp')
    
    install.packages(c('MASS', 'akima', 'robustbase', 'cobs', 'robust', 'mgcv', 'scatterplot3d', 'quantreg', 'rrcov', 'lars', 'pwr', 'trimcluster', 'mc2d', 'psych', 'Rfit', 'DepthProc', 'class', 'fda', 'e1071', 'rankFD'))
    devtools::install_github("jerryzhujian9/WRS_WRSC", subdir="pkg")
