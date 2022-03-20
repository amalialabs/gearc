FROM openjdk:17-alpine

RUN apk update
RUN apk upgrade
RUN apk add R=4.0.1

ADD secretlab1 secretlab1/

