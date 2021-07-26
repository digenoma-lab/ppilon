################# BASE IMAGE #####################
FROM nfcore/base

#FROM ubuntu:20.04
##site to test docker configuration files
# https://labs.play-with-docker.com/
################## METADATA #######################
LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="ppilon"
LABEL software.version="0.1"
LABEL about.summary="Container image containing all requirements for polishing assemblies"
LABEL about.home="https://github.com/camoragaq/"
LABEL about.documentation="https://github.com/camoragaq/ppilon/README.md"
LABEL about.license_file="https://github.com/camoragaq/ppilon/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER Carol Moraga <camoragaq@gmail.com>
################## INSTALLATION ######################
#polishing tools
COPY environment-ppilon.yml /
RUN conda env create -n ppilon -f /environment-ppilon.yml && conda clean -a
#ENV PATH /miniconda/envs/ppilon/bin:$PATH #env ubuntu
#ENV for nfcore-base
ENV PATH /opt/conda/envs/ppilon/bin:$PATH

