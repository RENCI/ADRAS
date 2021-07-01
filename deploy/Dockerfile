##############
# Docker file for the creation of the adcirc support tools.
#
# to create image: docker build -t adcirc_supp:latest .
# to push image:
#       docker tag adcirc_supp:latest renciorg/adcirc_supp:latest
#       docker push renciorg/adcirc_supp:latest
##############
FROM continuumio/miniconda3

# get some credit
LABEL maintainer="bblanton@renci.org"

# make sure the container is up-to-date
RUN apt-get --assume-yes update
RUN apt-get --assume-yes install vim bc

# Initial python install
RUN conda install  -y python=3.7

# update conda and set the download channel
RUN conda update conda && \
    conda config --add channels conda-forge

# tell conda what the shell is
RUN conda init bash

# make d directory for the repos and go there
RUN mkdir /repo
WORKDIR /repo

# get the repos
RUN git clone --branch hazus https://github.com/RENCI/ADRAS.git

# move to the install directory
WORKDIR ADRAS/deploy

# install the conda requirements
RUN conda install --yes --file conda_requirements.txt

# install the pip requirements
#RUN pip install -r pip_requirements.txt

# change to the noaa directory
#WORKDIR /repo/noaa_coops

# run the noaa seteup
#RUN python setup.py install --user

# change to the pipelines directory                        
WORKDIR /repo/ADRAS

# make a temporary directory for the output.
# this should eventually point a common PV
RUN mkdir work

# set the python path
ENV PYTHONPATH=/repo/ADRAS

# set the location of the output directory
ENV RUNTIMEDIR=/data

# lastly, get the aws credentials file
#RUN mkdir -p /repo/aws
#COPY aws_config /repo/aws/config
#COPY aws_credentials /repo/aws/credentials
#ENV AWS_CONFIG_FILE=/repo/aws/config
#ENV AWS_SHARED_CREDENTIALS_FILE=/repo/aws/credentials 


