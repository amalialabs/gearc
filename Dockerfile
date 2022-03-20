FROM openjdk:17-alpine

RUN apk update
RUN apk upgrade
RUN apk add R

ADD secretlab1 secretlab1/
ADD secretlab.jar secretlab.j

