FROM ubuntu:14.04
MAINTAINER Dan Lanier <mradmstr514226508@gmail.com>

# Docker image for Jupyterhub distribution of https://github.com/dlanier/FlyingMachineFractal.git

RUN apt-get update && apt-get install -y vim python3-pip git libblas-dev liblapack-dev libatlas-base-dev gfortran libfreetype6-dev libxft-dev python3-setuptools

RUN pip3 install -I numpy==1.11.1 pandas==0.18.1 scipy==0.19.1 scikit-learn==0.17.1 pyyaml xmlrunner

RUN pip3 install --upgrade pip

RUN pip3 install jupyter
RUN pip3 install ipywidgets

RUN jupyter nbextension enable --py --sys-prefix widgetsnbextension

RUN pip3 install jupyterhub
RUN apt-get update && apt-get install -y npm nodejs-legacy
RUN npm install -g configurable-http-proxy

RUN pip3 install jupyter_contrib-nbextensions
RUN jupyter contrib nbextension install --user

# RUN pip3 install --upgrade Pillow

# Set up working directory
ENV SRC_LOC /home
WORKDIR ${SRC_LOC}

# Oct 29 2017 time_to_crash at Pillow = 22m57.675s
#
# post build change:
#	move Pillow to end
#
# rebuild from intermediate container:
#	aka1spanky/fractal_jupyterhub:base_hub_10_29_2017

