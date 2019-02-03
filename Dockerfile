# Base image https://hub.docker.com/u/rocker/
FROM rocker/tidyverse:3.4.1

## Install extra R packages using requirements.R
## Specify requirements as R install commands e.g.
## 
## install.packages("<myfavouritepacakge>") or
## devtools::install("SymbolixAU/googleway")

# install packrat
# RUN R -e 'install.packages("packrat", repos="http://cran.rstudio.com", dependencies=TRUE, lib="/usr/local/lib/R/site-library");'
# RUN R -e 'install.packages("data.table", repos="http://cran.rstudio.com", dependencies=TRUE, lib="/usr/local/lib/R/site-library");'

# COPY ./DockerConfig/requirements.R /tmp/requirements.R 
# RUN Rscript /tmp/requirements.R

# # uncomment to include shiny server
# #RUN export ADD=shiny && bash /etc/cont-init.d/add

# create an R user
ENV USER rstudio

# pull current version of blm replication
# COPY --chown=rstudio:rstudio ./ /home/$USER/blm-replicate
# COPY ./ /home/$USER/blm-replicate

# run packrat - get all dependencies
WORKDIR /home/$USER
RUN git clone https://github.com/tlamadon/blm-replicate.git 
#RUN R -e 'packrat::restore("blm-replicate");'
RUN R -e 'devtools::install_deps("blm-replicate");'
RUN R -e 'devtools::install_github("setzler/textables");'

# build the blm-replicate library
RUN R CMD INSTALL --no-multiarch --with-keep.source blm-replicate

# make it visible to rstudio
RUN chown -R rstudio:rstudio .

