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
LONGOPTS=outdir:,genelist:,obo:,mapping:


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
obo=
mapping=
outdir=
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
      --)
        shift
        break
        ;;
      *)
        shift
        ;;
    esac
done

# handle non-option arguments
if [[ $# -ne 0 ]]; then
    echo "$0: empty flag detected, is this intentional?!"
    #exit 4
fi

mappingCall=
oboCall=

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
  outCall="--out $outdirPath"
fi

winpty docker run --pull=always $genelist $obo $mapping $outTmp --rm -it hadziahmetovic/secretlab1 secretlab \
  $genelistCall $oboCall $mappingCall $outCall


# ATTENTION!!! git bash must have enabled symlinks, look in c/programdata/git/config for [core] -> symlinks=true