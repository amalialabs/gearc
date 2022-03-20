FROM openjdk:17-alpine

RUN apk update
RUN apk upgrade
RUN apk add R

ADD secretlab1 secretlab1/
RUN mv /secretlab1/secretlab.j /secretlab1/secretlab.jar

RUN mkdir /out
ADD data data/

