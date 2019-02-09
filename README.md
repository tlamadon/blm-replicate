# Replication package for BLM

This git repository contains all the code to replicate the results of Bonhomme Lamadon and Manresa "A distributional Framework for matched employer-employee data". The working-paper version is available [here](http://lamadon.com/paper/blm.pdf). Virtually all code is based on the R platform. 

If you are looking for the R package to use the method, you should look at the [rblm package](https://tlamadon.github.io/rblm/index.html). It includes most of the estimators available here, and we keep updating it.

The present replication package is built as an R package that can be easily installed on any system. All package dependencies are handled using [packrat](https://rstudio.github.io/packrat/). This guarantees that results can be reproduced using the exact versions of all the libraries that were used at the time the paper was written. We also provide a Docker container to ensure full portability. 

Importantly, reproducing the results on Swedish data __requires access to the administrative data from Sweden__. Researchers need to apply to get access to such data. We recommend contacting the [IFAU](https://www.ifau.se/). The institute is hosting this replication package that can be accessed and ran on the data on their servers. The reference name for our project is `IFAU-2015-65`.



## How do I run this?

The simplest way to use this replication package is to rely on the docker container that we have created. 

### Solution 1: get it running in less than 10 minutes, run our docker container

If you want to see the code in action, you only need the following steps:

1.  install the [docker app](https://www.docker.com/get-started) (if you don't have it already)
2. download and start our replication container with only one command: `docker run -d --rm -e PASSWORD=blm -p 8787:8787 tlamadon/blm-replicate`
3. open your browser at [http://localhost:8787](http://localhost:8787/) and use login `rstudio` and password `blm` 
4. finally run `source("inst/main.R")`

This will give you access to a fully functioning RStudio with the installed libraries and the code necessary to run the replication code. By default, this will run all of the code using a __synthetic data set__. See later how to get access to Swedish data, and load it into the container.

Note: make sure the docker app does not limit memory access to less than 16Gb. See [here](https://stackoverflow.com/questions/44417159/docker-process-killed-with-cryptic-killed-message). 

### Solution 2: install the replication package into your R environment

If you have your own running R system, and you want to run this replication package in your environment, you can directly install the package. In this case we recommend that you make use of the [packrat](https://rstudio.github.io/packrat/) configuration we are providing.

1. Download the replication package. We recommend to simply clone the github repository, ie:  `git clone https://github.com/tlamadon/blm-replicate.git`
2. Start R inside the replication package.

In R, run the following commands:

```R
# installing the package locally in your R env.

install.packages("pakcrat") # make sure that packrat is available
install.packages("devtools") # make sure that devtools is available

source("packrat/init.R") # initialize the packrat environment
packrat::restore()       # make sure all is up to date

devtools::install(".")   # build the replication package

source("inst/main.R")    # fire up the replication
```

## Overview of the replication package

The main enty point is [inst/main.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/main.R). It will __automatically__ run all the necessary steps in the other files in order to reproduce all the results of the paper. Note however that this would take a very long time as it will starts some bootstraps procedure. The code will generate all figures and tables and put them into folder called `tmp`  by default.

We invite resesearchers to read through [inst/main.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/main.R) which has explicit calls for each subsets of the paper. 

### Organization of the code

 - All the heavy lifting such as the estimators and simulation codes are in the [R/*.r](https://github.com/tlamadon/blm-replicate/tree/master/R) folder. This is the usual way to store functions in an R package.
 - [inst/server/data-section-static.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/data-selection-static.r) contains the code that __processes the data inputs__ to prepare the data for the static estimation.
 - [inst/server/data-section-dynamic.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/data-selection-dynamic.r) contains the code that __processes the data inputs__ to prepare the data for the dynamic estimation.
 - [inst/server/estimation-static.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/estimation-static.r) contains the code that runs the estimations for the static version of the model 
 - [inst/server/estimation-dynamic.r](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/estimation-dynamic.r) contains code that runs the different __estimations__ for the dynamic version of the model.
 - [inst/server/fig-blm.R](https://github.com/tlamadon/blm-replicate/blob/master/inst/server/fig-blm.R) contrains functions that generate all of the __figures and tables__ in the paper.
 - Finally we are collecting the sequence of calls that generates all the results 

## Replicating the results on Swedish data

### 

The main data source should be the following list of:

- `selectedf0educ1.dta` 
- `selectedf0educ2.dta` 
- `selectedf0educ3.dta` 
- `selectedf1educ1.dta` 
- `selectedf1educ2.dta` 
- `selectedf1educ3.dta` 
- `selectedfirms9708.dta`

and should be available to the project as a relative path `../data/`. Researchers can apply for access at the [IFAU](https://www.ifau.se/).

 




