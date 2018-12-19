# Replication package for BLM

This git repository contains all the code to replicate the results of Bonhomme Lamadon and Manresa "A distributional Framework for matched employer-employee data". The working-paper version is available [here](http://lamadon.com/paper/blm.pdf). Virtually all code is based on the R platform. 

If you are looking for the R package to use the method, you should look at the [rblm package](https://tlamadon.github.io/rblm/index.html). It includes most of the estimators available here, and we keep updating it, whereas we keep the code here to the version we used to generate the paper.

The present replication package is build as an R package that can be easily installed on any system. All package dependencies are handled using [packrat](https://rstudio.github.io/packrat/). This guarantees that results can be reproduced using the exact versions of all the libraries that were used at the time the paper was written.

Importantly, reproducing the results on Swedish data requires access to the administrative data from Sweden. Researcher have to apply to get access to such data. We recommend contacting the [IFAU](https://www.ifau.se/). The institute is hosting this replication package that can be accessed and ran on the data on their servers with their data.

# Overview of the replication package

## The data source: the swedish registry

The main data source should be the following list of:

- `selectedf0educ1.dta` 
- `selectedf0educ2.dta` 
- `selectedf0educ3.dta` 
- `selectedf1educ1.dta` 
- `selectedf1educ2.dta` 
- `selectedf1educ3.dta` 
- `selectedfirms9708.dta`

and should be available to the project as a relative path `../data/`. Researchers can appy for access at the IFAU.

## The main replication file

The main replication file is [inst/main.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/main.R).  It will run all the necessary steps in order to reproduce all the results of the paper. Note however that running through it would take several weeks even when using 15 cpus. 




