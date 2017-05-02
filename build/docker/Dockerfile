FROM ubuntu:14.04
MAINTAINER D Lanier <lanier4@illinois.edu>

RUN apt-get update && \
	apt-get install -y vim python3-pip git libblas-dev liblapack-dev libatlas-base-dev gfortran libfreetype6-dev libxft-dev

RUN pip3 install -I numpy==1.11.1 pandas==1.18.1 scipy==0.18.0 scikit-learn==0.17.1 pyyaml