FROM openjdk:17-alpine

RUN apk update
RUN apk upgrade
RUN apk add R

ADD secretlab1 secretlab1/
RUN mv secretlab.jar secretlab1/secretlab.j

RUN mkdir /out
ADD data data/

