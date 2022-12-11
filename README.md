# *gEArc* v1.0

# Getting started
For convenience, the *gEArc* tool is provided as an image. <br>
The latest version is available under:<br>
- [amalialabs/gearc](https://github.com/amalialabs/gearc)


## Description
gEArc is a gene set enrichment tool correcting for several flaws of standard over-representation analysis with Gene Ontology terms. It operates only on a clear measured subset of the input gene list and avoids the classic binary threshold for significance while providing robust enrichment results. Instead of arbitrary fold-change and FDR cutoffs, gEArc extends them to an allowed interval so that edge cases of measured genes are evaded. Based on those intervals, genes are weighted and assigned into three categories: SIG-CORE, FLEX-CORE, SIGNON-CORE. From these categories, genes are then sampled several times for the robust enrichment. This way, the FLEX-CORE includes genes for which no clear statement about their relevance for the enrichment can be made. In the end, the robust enrichment provides the user with significant GOO terms enriched at a default 95% quantile.


## Usage

### Starting gEArc
The easiest way to start is by executing either `run.sh` or `run_windows.sh`
depending on whether you are using unix or windows. Both scripts make it easier
to map all of the needed variables to the container in which the computation will
take place.

	run.sh --genelist <file> [--obo <file>] [--mapping <file>] [--n <integer>] [--root <GO root>] [--FDR <float>] [--FC <float>] [--expectedChange <LOW|AVERAGE|HIGH>]

    --genelist <file>
      The genelist should be in the format of gene_id, log2 fold change, and fdr.
      If you have a different format you may use the reformatGeneList.sh script in
      order to get the desired input format.

    [--obo <file>]
      `default go-basic.obo from Gene ontology`
	  obo file containing the GO DAG structure. Use only if you want to use a customized ontology.
      Make sure the format corresponds to the obo specification as used by Gene Ontology.

    [--mapping <file>]
      `default ensembl GO gene mapping release 106 - human`
	  mapping file containing Gene to GO Node mappings. Use only if you want to use a customized gene network.
      The mapping extracted from Ensembl and external references to GO. File format is a tsv with the columns:
      ensembl_id, hgnc, and pipe (|) separated GO_ids.


    [--n <integer>]
      `default 50`
	  set the number of iterations for the subsampling.


    [--root <GO root>]
      `default biological_process`
	  select the GO root for which the analysis should be performed. Includes biological_process,
      molecular_function, and cellular_component.


    [--FDR <float>]
      `default 0.01`
	  select FDR cutoff for the classification of significance.


    [--FC <float>]
	  `default 1.0 (lfc)`
      select FC cutoff for the classification of significance.

    [--expectedChange <LOW|AVERAGE|HIGH>]
      `default: AVERAGE`.
	  Select the parameter based on the expected differences in your data. Viral infections or stress
      induction generally result in HIGH changes. AVERAGE changes are to be expected from time series or general changes
      between different tissues. LOW changes are expected from minor alterations such as specific knock outs (KO) where only
      specific networks are affected.

### Requirement ###
#### Prepare required input data: <br>
- differential gene expression results file in the following tsv-format <br>

Gene   |  FDR    |  LFC
-------------|-------|-----
ENSG00000105     |  0.15   |  1.1
ENSG00000105939  |  0.7    |  0.2

#### Docker or Podman
The `run.sh` script will attempt to pull the image from our dockerhub repository.
For this you will need to have either podman or docker installed, default is podman, but you may change to docker in the `run.sh` on the last line.

The `run_windows.sh` script uses only docker.


## Example
### Unix
$ run.sh --genelist DESeq_hisat.tsv --outdir ./out --expectedChange HIGH --n 100 > run.log

### Windows
Important: this is not for the linux subsystem, but the git shell.
$ run_windows.sh --genelist DESeq_hisat.tsv --outdir ./out --expectedChange HIGH --n 100 > run.log



## Files

### Genelist format conversion
The genelist should be in the format of Gene_id, LFC (log2 fold change), and FDR.
If you have a different format you may use the `reformatGeneList.sh` script in order to get the desired input format.
It takes 4 mandatory inputs: `--gene` `--lfc` `--fdr` `--input`.
If the `--out` param is set then the output is written to that file, else
it is returned to standard out.
```shell script
$ head DESeq_hisat.out 
ENE.ID log2FC  RAW.PVAL        ADJ.PVAL
ENSG00000237973 -3.89663231733891       6.86535048511563e-05    0.000294971690634794
ENSG00000228794 -0.700529152285993      0.455927017552331       0.562476661534651
ENSG00000223764 -1.59236001133306       0.0454846640889237      0.0884531248524819
ENSG00000187634 0.966894090405194       0.236147768382723       0.335949929405359

$ ./reformatGeneList.sh --gene 1 --lfc 2 --fdr 4 --input DESeq_hisat.out  --out DESeq_hisat.out .reformatted

$ head ../data/corona.reformatted 
gene    diffexp.log2fc  diffexp.fdr
ENSG00000143067 -2.230  0.000e+00
ENSG00000168394 -3.130  0.000e+00
ENSG00000181381 -3.050  0.000e+00
ENSG00000204592 -1.390  0.000e+00
```

### Mapping
wget http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
gunzip Homo_sapiens.GRCh38.106.gtf.gz
cat Homo_sapiens.GRCh38.103.gtf | grep -P "\tgene\t" | sed 's/^.*gene_id "//' | sed 's/"; gene_.*$//' | sort -u > Homo_sapiens.GRCh38.103.gtf.geneIDs
