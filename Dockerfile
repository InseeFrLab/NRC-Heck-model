FROM rocker/verse:3.6.3

RUN apt-get update \
    && apt upgrade -y

RUN Rscript -e 'update.packages()'
RUN Rscript -e 'install.packages("sampleSelection", version = "1.2-12")'
RUN Rscript -e 'install.packages("latex2exp", version = "0.5.0")'
