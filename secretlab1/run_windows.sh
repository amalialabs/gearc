#!/bin/bash -x

echo $@
params=("$@")

# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset

# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I’m sorry, `getopt --test` failed in this environment.'
    exit 1
fi

OPTIONS=
LONGOPTS=outdir:,genelist:,obo:,mapping:,root:,expectedChange:,n:,FDR:,FC:

usage () {
    cat <<HELP_USAGE
$0  --genelist <file> [--obo <file>] [--mapping <file>] [--n <integer>] [--root <GO root>] [--FDR <float>] [--FC <float>] [--expectedChange <LOW|AVERAGE|HIGH>]

    --genelist <file>
      The genelist should be in the format of gene_id, log2 fold change, and fdr.
      If you have a different format you may use the reformatGeneList.sh script in
      order to get the desired input format.

    [--obo <file>]
      obo file containing the GO DAG structure. Use only if you want to use a customized ontology.
      Make sure the format corresponds to the obo specification as used by Gene Ontology.

    [--mapping <file>]
      mapping file containing Gene to GO Node mappings. Use only if you want to use a customized gene network.
      The mapping extracted from Ensembl and external references to GO. File format is a tsv with the columns:
      ensembl_id, hgnc, and pipe (|) separated GO_ids.


    [--n <integer>]
      set the number of iterations for the subsampling.


    [--root <GO root>]
      select the GO root for which the analysis should be performed. Includes biological_process,
      molecular_function, and cellular_component.


    [--FDR <float>]
      select FDR cutoff for the classification of significance.


    [--FC <float>]
      select FC cutoff for the classification of significance.

    [--expectedChange <LOW|AVERAGE|HIGH>]
      default: AVERAGE. Select the parameter based on the expected differences in your data. Viral infections or stress
      induction generally result in HIGH changes. AVERAGE changes are to be expected from time series or general changes
      between different tissues. LOW changes are expected from minor alterations such as specific knock outs (KO) where only
      specific networks are affected.
HELP_USAGE
}

# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"
genelist=
obo=
mapping=
outdir=
n=
FDR=
FC=
root=
expectedChange=
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
		  --outdir)
        outdir="$2"
        shift 2
        ;;
		  --genelist)
        genelist="$2"
        shift 2
        ;;
      --obo)
        obo="$2"
        shift 2
        ;;
      --mapping)
        mapping="$2"
        shift 2
        ;;
      --n)
        n="$2"
        shift 2
        ;;
      --root)
        root="$2"
        shift 2
        ;;
      --FDR)
        FDR="$2"
        shift 2
        ;;
      --FC)
        FC="$2"
        shift 2
        ;;
      --expectedChange)
        expectedChange="$2"
        shift 2
        ;;
      --)
        shift
        break
        ;;
      *)
        shift
        ;;
    esac
done

if [[ "$genelist" == "" ]] ; then
    echo STDERR "no genelist provided. Program will exit."
    usage
    exit 0
fi

# handle non-option arguments
if [[ $# -ne 0 ]]; then
    echo "$0: empty flag detected, is this intentional?!"
    usage
    #exit 4
fi

mappingCall=
oboCall=
outCall=
nCall=
outTmp=

if [[ "$mapping" != "" ]]; then
  mappingPath=$(echo $(readlink -f $mapping) | sed 's/\/\\/g')
  mapping="-v $mappingPath:/input/mapping.tsv";
  mappingCall="--mapping //input/mapping.tsv"
fi
if [[ "$obo" != "" ]]; then
  oboPath=$(echo $(readlink -f $obo) | sed 's/\/\\/g')
  obo="-v $oboPath:/input/go.obo";
  oboCall="--obo //input/go.obo"
fi

genelistPath=$(echo $(readlink -f $genelist) | sed 's/"\"/"\\"/g')
genelist="-v /$genelistPath:/input/genelist.tsv"
genelistCall="--genelist //input/genelist.tsv"

if [[ "$outdir" != "" ]]; then
  mkdir -p $outdir
  outdirPath=$(echo $(readlink -f $outdir) | sed 's/"\"/"\\"/g')
  outTmp="-v /$outdirPath:/$outdirPath"
  outCall="--out /$outdirPath"
fi

if [[ "$n" != "" ]]; then
  nCall="--n $n"
fi

rootCall=
FDRCall=
FCCall=
expectedChangeCall=

if [[ "$root" != "" ]]; then
  rootCall="--root $root"
fi

if [[ "$FDR" != "" ]]; then
  FDRCall="--FDR $FDR"
fi

if [[ "$FC" != "" ]]; then
  FCCall="--FC $FC"
fi

if [[ "$expectedChange" != "" ]]; then
  expectedChangeCall="--expectedChange $expectedChange"
fi


winpty docker run --pull=always $genelist $obo $mapping $outTmp --rm -it amalialabs/gearc secretlab \
  $genelistCall $oboCall $mappingCall $outCall $nCall $rootCall $FDRCall $FCCall $expectedChangeCall


# ATTENTION!!! git bash must have enabled symlinks, look in c/programdata/git/config for [core] -> symlinks=true