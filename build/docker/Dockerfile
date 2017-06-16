# https://github.com/jupyter/docker-stacks
FROM jupyter/scipy-notebook
MAINTAINER Dan Lanier <lanier4@illinois.edu>

USER root

# Install packages and their dependencies
RUN apt-get update && apt-get install -y python3-pip libblas-dev liblapack-dev libatlas-base-dev gfortran libfreetype6-dev libxft-dev
RUN apt-get update && apt-get install python3-setuptools

RUN conda install numpy==1.11.1 pandas==0.18.1 pyyaml Pillow IPython

USER $NB_USER

