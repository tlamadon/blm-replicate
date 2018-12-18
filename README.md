# Replication package for BLM


This git repository contains all the code to replicate the results of Bonhomme Lamadon and Manresa "A distributional Framework for matched data". Virtually all code is based on the R platform. 

The package is build as an R package that can be easily installed on any R system. All package dependencies are handled using `packrat`. This guarantees that results can be reproduced using the exact versions that were used at the time the paper was written.

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



## The main replication file

The main replication file is `inst/main.r`.  It will run all the necessary steps in order to reproduce all the results of the paper. Note however that running through it would take several weeks even when using 15 cpus. 

