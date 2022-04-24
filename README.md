# *secretlab1* v1.0

# Getting started
For convenience, the *secretlab1* tool is provided as an image. <br>
The latest version is available under:<br>
- [secretlab/project1](INSERT link to repo here)

### Step1 ###
Prepare needed input data: <br>
- differential results file in the following tsv-format <br>

Gene   |  fdr    |  lfc
-------------|-------|-----
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

## Files

#### Genelist
The genelist should be in the format of gene_id, log2 fold change, and fdr.
If you have a different format you may use the `reformatGeneList.sh` script in order to get the desired input format.
It takes 4 mandatory inputs: `--gene` `--lfc` `--fdr` `--input`.
If the `--out` param is set then the output is written to that file, else
it is returned to standard out.
```shell script
$ head ../data/corona.empires.outECC 
gene    diffexp.fdr     diffexp.log2fc  diffsplic.most.signif.test      diffsplic.fdr   diffsplic.difflog2fc
ENSG00000143067 0.000e+00       -2.230  untested        1.000   0.000e+00
ENSG00000168394 0.000e+00       -3.130  ENSG00000168394.merged.ENST00000354258.ENST00000643049_VS_excl.ENST00000643049  0.240   -0.919
ENSG00000181381 0.000e+00       -3.050  ENSG00000181381.merged.ENST00000510590.ENST00000511577_VS_excl.ENST00000510590  0.265   0.862
ENSG00000204592 0.000e+00       -1.390  ENSG00000204592.excl.ENST00000484194_VS_excl.ENST00000493699    0.060   -1.013

$ ./reformatGeneList.sh --gene 1 --lfc 3 --fdr 2 --input ../data/corona.empires.outECC --out ../data/corona.reformatted

$ head ../data/corona.reformatted 
gene    diffexp.log2fc  diffexp.fdr
ENSG00000143067 -2.230  0.000e+00
ENSG00000168394 -3.130  0.000e+00
ENSG00000181381 -3.050  0.000e+00
ENSG00000204592 -1.390  0.000e+00
```

## Starting secretlab
The easiest way to start is by executing either `run.sh` or `run_windows.sh`
depending on whether you are using unix or windows. Both scripts make it easier
to map all of the needed variables to the container in which the computation will
take place.



### Mapping
wget http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
gunzip Homo_sapiens.GRCh38.106.gtf.gz
cat Homo_sapiens.GRCh38.103.gtf | grep -P "\tgene\t" | sed 's/^.*gene_id "//' | sed 's/"; gene_.*$//' | sort -u > Homo_sapiens.GRCh38.103.gtf.geneIDs