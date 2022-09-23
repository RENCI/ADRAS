##############
# Docker file for the creation of the ADRAS/HAZUS geotiff files.
# stage 1: get a conda virtual environment
##############
FROM continuumio/miniconda3 as build

# get some credit
LABEL maintainer="Brian Blanton (bblanton@renci.org)"

# extra metadata
LABEL version="v0.0.1"
LABEL description="ADRAS/HAZUS Dockerfile."

# update conda
RUN conda update conda

# Create the virtual environment
COPY environment.yml .
RUN conda env create -f environment.yml

# install conda pack to compress this stage
RUN conda install -c conda-forge conda-pack

# conpress the virtual environment
RUN conda-pack -n adras-hazus -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# fix up the paths
RUN /venv/bin/conda-unpack

##############
# stage 2: create a python implementation using the stage 1 virtual environment
##############
FROM python:3.9-slim

RUN apt-get update

# install wget and bc
RUN apt-get install -y wget bc

# clear out the apt cache
RUN apt-get clean

# add user nru and switch to it
RUN useradd --create-home -u 1000 nru
USER nru

# Create a directory for the log
RUN mkdir -p /home/nru/ADRAS/logs

# move the the code location
WORKDIR /home/nru/ADRAS

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# make the virtual environment active
ENV VIRTUAL_ENV /venv
ENV PATH /venv/bin:$PATH

# Copy in the rest of the code
COPY config config
COPY misc misc
COPY SRC SRC
COPY utilities utilities
COPY ./*.sh ./

# set the python path
ENV PYTHONPATH=/home/nru/ADRAS

# set the location of the output directory
ENV RUNTIMEDIR=/data
ENV PKLDIR=/data/pkldir

# set the log dir. use this for debugging if desired
ENV LOG_PATH=/data/logs

# example command line
# bash compute_geotiffs.sh http://tds.renci.org:80/thredds/fileServer/2022/al07/24/NCSC_SAB_v1.23/hatteras.renci.org/ncsc123-al07-sb55.01/nhcOfcl
