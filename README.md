# Replication package for BLM

This git repository contains all the code to replicate the results of Bonhomme Lamadon and Manresa "A distributional Framework for matched employer-employee data". The working-paper version is available [here](http://lamadon.com/paper/blm.pdf). Virtually all code is based on the R platform. 

If you are looking for the R package to use the method, you should look at the [rblm package](https://tlamadon.github.io/rblm/index.html). It includes most of the estimators available here, and we keep updating it, whereas we keep the code here to the version we used to generate the paper.

The present replication package is build as an R package that can be easily installed on any system. All package dependencies are handled using [packrat](https://rstudio.github.io/packrat/). This guarantees that results can be reproduced using the exact versions of all the libraries that were used at the time the paper was written.

Importantly, reproducing the results on Swedish data requires access to the administrative data from Sweden. Researcher have to apply to get access to such data. We recommend contacting the [IFAU](https://www.ifau.se/). The institute is hosting this replication package that can be accessed and ran on the data on their servers with their data.

# Overview of the replication package

### Organization of the code

 - All the heavy lifting such as the estimators and simulation codes are in the [R/*.r](https://github.com/tlamadon/blm-replicate/tree/master/R) folder. This is the usual way to store functions in an R package.
 - [inst/server/data-section-static.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/data-selection-static.r) contains the code that __process the data inputs__ to prepare the data for the static estimation.
 - [inst/server/data-section-dynamic.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/data-selection-dynamic.r) contains the code that __process the data inputs__ to prepare the data for the dynamic estimation.
 - [inst/server/estimation-static.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/estimation-static.r) contains the code that runs the __estimations__ for the static version of the model 
 - [inst/server/estimation-dynamic.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/estimation-dynamic.r) contains code tha run the different __estimations__ for the dynamic version of the model.
 - [inst/server/fig-blm.R](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/fig-blm.R) contrains functions that generate all of the __figures and tables__ in the paper.
 - Finally we are collecting the sequence of calls that generates all the results in [inst/main.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/main.R). It will __run all the necessary steps__ in order to reproduce all the results of the paper. 

### The data source: the swedish registry

The main data source should be the following list of:

- `selectedf0educ1.dta` 
- `selectedf0educ2.dta` 
- `selectedf0educ3.dta` 
- `selectedf1educ1.dta` 
- `selectedf1educ2.dta` 
- `selectedf1educ3.dta` 
- `selectedfirms9708.dta`

and should be available to the project as a relative path `../data/`. Researchers can appy for access at the IFAU.

 




