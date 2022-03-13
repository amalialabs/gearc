# *secretlab1* v1.0

# Getting started
For convenience, the *secretlab1* tool is provided as an image. <br>
The latest version is available under:<br>
- [secretlab/project1](INSERT link to repo here)

### Step1 ###
Prepare needed input data: <br>
- differential results file in the following format <br>
Gene             |  fdr    |  lfc
-------------    | ------- | -----
ENSG00000105     |  0.15   |  1.1
ENSG00000105939  |  0.7    |  0.2

### Step 2 ###
If you have podman or docker ready to run, then simply execute: <br>
```shell script
podman run -v diffExp.tsv:/input/diffExp.tsv secretlab/project1 [--plots] [--ensembl|gaf]
```

### Step 3 / Final step
*Enjoy life because you have great results :)*

### Build it yourself
To create the image simply clone the repository and from within run: <br>
```shell script
podman build -t {preffered_name} .
```
where the `-t` specifies the tag the image will have and the `.` 
specifies the *context*, i.e. where the *Dockerfile* is located.
This will perform all instructions from the *Dockerfile* and thus assemble
the runnable image.

#### Dockerfile
Typically a *Dockerfile* starts of with a *base-image*, in our case
this is preferrably a Java or R base.
```shell script
FROM openjdk:17-alpine

RUN apk update
RUN apk upgrade
RUN apk add R=4.0.1

ADD scripts /home/scripts
```

or alternatively

```shell script
FROM rocker/r-ver:4.0.1

RUN apt-get update && \
    apt-get install -y openjdk-17-jdk && \
    apt-get install -y liblzma-dev && \
    apt-get install -y libbz2-dev

RUN Rscript -e "install.packages('rJava')"
```