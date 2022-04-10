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

## for podman users replace here docker with 'podman'
## for windows users replace here docker with 'winpty docker'
docker run --pull=always $obo $mapping -v $outdir:/out/ -v $genelistPath:$genelistPath --rm -it hadziahmetovic/secretlab1 secretlab \
  $genelistCall $oboCall $mappingCall