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
LONGOPTS=gene:,fdr:,lfc:,input:,out:


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

out=
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
		  --gene)
        gene="$2"
        shift 2
        ;;
		  --fdr)
        fdr="$2"
        shift 2
        ;;
      --lfc)
        lfc="$2"
        shift 2
        ;;
      --input)
        input="$2"
        shift 2
        ;;
      --out)
        out="$2"
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

if [[ "$out" != "" ]]; then
  cat $input | awk -v gene=$gene -v fdr=$fdr -v lfc=$lfc '{print $gene "\t" $lfc "\t" $fdr}' > $out
else
  cat $input | awk -v gene=$gene -v fdr=$fdr -v lfc=$lfc '{print $gene "\t" $lfc "\t" $fdr}'
fi
