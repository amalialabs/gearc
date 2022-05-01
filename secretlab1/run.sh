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
    #exit 4
fi

mappingCall=
oboCall=
nCall=
outdirCall=

if [[ "$mapping" != "" ]]; then
  mappingPath=`readlink -f $mapping`;
  mapping="-v $mappingPath:$mappingPath";
  mappingCall="--mapping $mappingPath"
fi
if [[ "$obo" != "" ]]; then
  oboPath=`readlink -f $obo`;
  obo="-v $oboPath:$oboPath";
  oboCall="--obo $oboPath"
fi

genelistPath=`readlink -f $genelist`
genelistCall="--genelist $genelistPath"

if [[ "$outdir" != "" ]]; then
  mkdir -p $outdir
  outdirPath=`readlink -f $outdir`
  outdirCall="--out $outdirPath"
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

## for docker users replace here podman with 'docker'
## for windows users replace here podman with 'winpty docker'
podman run --pull=always $obo $mapping -v $outdir:/out/ -v $genelistPath:$genelistPath --rm -it hadziahmetovic/secretlab1 secretlab \
  $genelistCall $oboCall $mappingCall $outdirCall $nCall $rootCall $FDRCall $FCCall $expectedChangeCall