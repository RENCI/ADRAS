##############
# Docker file for the creation of the adcirc support tools.
#
# to create image: docker build -t adcirc_supp:latest .
# to push image:
#       docker tag adcirc_supp:latest renciorg/adcirc_supp:latest
#       docker push renciorg/adcirc_supp:latest
##############
FROM continuumio/miniconda3

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

RUN useradd -m -d /home/nru -u 1000 nru
RUN chmod 777 /home/nru
USER nru

# make d directory for the repos and go there
RUN mkdir /home/nru/repo
WORKDIR /home/nru/repo

# get the repos
#RUN git clone --branch rotation https://github.com/RENCI/ADRAS.git
RUN git clone https://github.com/RENCI/ADRAS.git

# move to the install directory
WORKDIR /home/nru/ADRAS/deploy
#WORKDIR /repo/ADRAS
# RUN ls
# install the conda requirements

RUN conda install --yes --file conda_requirements.txt

# install the pip requirements
#RUN pip install -r pip_requirements.txt

# change to the noaa directory
#WORKDIR /repo/noaa_coops

# run the noaa seteup
#RUN python setup.py install --user

WORKDIR /home/nru/repo/ADRAS

# set the python path
ENV PYTHONPATH=/home/nru/repo/ADRAS

# set the location of the output directory
ENV RUNTIMEDIR=/data
ENV PKLDIR=/data/pkldir

