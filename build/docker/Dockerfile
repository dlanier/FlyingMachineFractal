FROM ubuntu:14.04
MAINTAINER Dan Lanier <mradmstr514226508@gmail.com>

# FlyingMachineFractal Docker BASE IMAGE does not include:   https://github.com/dlanier/FlyingMachineFractal.git

RUN apt-get update && apt-get install -y vim python3-pip git libblas-dev liblapack-dev libatlas-base-dev gfortran libfreetype6-dev libxft-dev python3-setuptools

#  apt-get install python-PIL
RUN apt-get install python-PIL

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

RUN pip3 install --upgrade Pillow
RUN pip3 install --upgrade matplotlib

# Set up working directory
ENV SRC_LOC /home
WORKDIR ${SRC_LOC}

