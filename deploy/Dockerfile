##############
# Docker file for the creation of the adcirc support tools.
#
# to create image: docker build -t adras:latest .
# to push image:
#       docker tag adras:latest renciorg/adras:latest
#       docker push renciorg/adras:latest
##############
FROM continuumio/miniconda3

# get some credit
LABEL maintainer="jtilson@renci.org"

# make sure the container is up-to-date
RUN apt-get update

# Start off getting Python 3.7.3
RUN conda create -name adras -y python=3.7.3 gdal rasterio pyyaml pandas numpy matplotlib netCDF4 geopandas

# tell conda what the shell is
RUN conda init bash

# make d directory for the repos and go there
RUN mkdir /repo
WORKDIR /repo

# get the repos
RUN git clone https://github.com/RENCI/ADRAS.git
               
# move to the install directory
WORKDIR ADRAS/INSTALL

# install the conda requirements
RUN conda install --yes --file conda_requirements.txt
                                
# install the pip requirements
RUN pip install -r pip_requirements.txt
                        
# change to the pipelines directory                        
WORKDIR /repo/ADRAS
                                   
# make a temporary directory for the output.
# this should eventually point a common PV
RUN mkdir temp

# set the python path
ENV PYTHONPATH=/repo/ADRAS

# set the location of the output directory
ENV RUNTIMEDIR=/data

