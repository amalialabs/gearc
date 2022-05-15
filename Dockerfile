FROM openjdk:17-alpine

RUN apk update
RUN apk upgrade
RUN apk add R
RUN apk add gcc g++ make
RUN apk add musl-dev
RUN apk add R-dev
RUN Rscript -e "install.packages('ggplot2', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('reshape2', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('tidyr', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('VennDiagram', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('data.table', repos='http://cran.rstudio.com/')"
#reshape2, VennDiagram, tidyr, data.table

RUN apk add maven
RUN apk add bash


ADD secretlab1 secretlab1/

RUN mkdir /out
ADD data data/

RUN cd secretlab1 && mvn clean package
RUN cd /secretlab1/target && mkdir libs && for file in `find /root/.m2/repository/ -name "*.jar"`; do cp $file libs/; done

RUN echo -e '#!/bin/bash \n java -cp "/secretlab1/target/libs/*":/secretlab1/target/ Handler "$@"' > /usr/bin/secretlab && \
    chmod +x /usr/bin/secretlab

