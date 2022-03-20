FROM openjdk:17-alpine

RUN apk update
RUN apk upgrade
RUN apk add R

ADD secretlab1 secretlab1/
RUN mv /secretlab1/secretlab.j /secretlab1/secretlab.jar

RUN mkdir /out
ADD data data/

RUN apk add maven
RUN cd secretlab1 && mvn package
RUN cd /secretlab1/target && mkdir libs && for file in `find /root/.m2/repository/ -name "*.jar"`; do cp $file libs/; done
RUN alias secretlab1='java -cp "/secretlab1/target/libs/*":`/secretlab1/target` Handler'

RUN Rscript -e "install.packages('ggplot2', dependencies=TRUE, repos='http://cran.rstudio.com/')"

