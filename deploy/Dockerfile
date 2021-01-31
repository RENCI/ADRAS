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
RUN conda install  -y python=3.7.3

# update conda and set the download channel
RUN conda update conda && \
    conda config --add channels conda-forge


# tell conda what the shell is
RUN conda init bash

# make d directory for the repos and go there
RUN mkdir /repo
WORKDIR /repo

# get the repos: Must be public or proviude user/pass
# Until the repo is public you must supply your github username:password
# 
RUN git clone https://username:password@github.com/RENCI/ADRAS.git
               
# move to the install directory
WORKDIR ADRAS/deploy

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

