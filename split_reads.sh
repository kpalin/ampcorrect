#!/bin/bash

#  Thu Jul 30 10:27:01 2015
#  kpalin@merit.ltdk.helsinki.fi

# For debugging
#set -o verbose 

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail





usage()  {
    echo -e "usage:\n$0 -F splitters.fasta sequence.(fast5|fasta)

Map short sequences in splitters.fasta to reads in sequence.fasta and
generate SAM like output file usable for e.g. split_fasta_by_sam.py

-F splitters.fast(a|5)  Split the reads according to alignment to these.
-d          Use default lastal arguments, instead of ones for nanopore reads
-h          Show this message and exit." >&2
    exit 1
    
}

#PATH=/mnt/cg8/projects/Nanopore/software/tmp/last-756/src:$PATH

function toFasta() {
    SEQFILENAME=$1
    FASTA=$2
    test ! -s ${FASTA}  # Check that output file doesn't exist


    isHDF=$(file ${SEQFILENAME}|grep -cF 'Hierarchical Data Format' || true)
    
    if [ ${isHDF} -eq 1 ]; then
	for READTYPE in 2D fwd rev
	do
	    poretools fasta --type ${READTYPE} ${SEQFILENAME} >${FASTA}
	    test ! -s ${FASTA} || break
	done
    else
	cp ${SEQFILENAME} ${FASTA}
    fi
}



while getopts  "F:hdS" flag
do
    case "$flag" in 
	F)
	    INPUT2=$(readlink -f "$OPTARG")
	    test -r ${INPUT2}
	;;
	S)
	    KEEPSAM=1
	    ;;
	d)
	    DEFAULTLASTAL=1
	    ;;
	h|*)
	usage
	;;
    esac
done
shift $((OPTIND-1)); OPTIND=1






trap usage EXIT

FAST5=$1


# Make temporary directory which will be cleaned at script exit
TEMPDIR=$(mktemp --directory )
function _cleanup {
    rm -r $TEMPDIR
}
trap _cleanup EXIT



BASE=$(basename ${FAST5} .fast5)

OUT=$(readlink -f ${BASE}.png)


FASTA1=${TEMPDIR}/${BASE}.fasta
MAF=${TEMPDIR}/${BASE}.maf


toFasta ${FAST5} ${FASTA1}


if [ -z ${INPUT2:-} ];
then
    FASTA2=${FASTA1}
else
    FASTA2=${TEMPDIR}/$(basename ${INPUT2} )_2.fasta
    toFasta ${INPUT2} ${FASTA2}

fi



SEQ1name=$(head -1 ${FASTA1}|sed -e 's/[>]\([^        ]\+\).*$/\1/')
SEQ2name=$(head -1 ${FASTA2}|sed -e 's/[>]\([^        ]\+\).*$/\1/')

SAM=$(readlink -f ${SEQ1name}_vs_${SEQ2name}.sam)


cd $TEMPDIR
lastdb  ${BASE} ${FASTA1}

if [ "${DEFAULTLASTAL:-0}" -eq 1 ];
then
    echo NO HERE
    lastal -T 0  ${BASE} ${FASTA2} >${MAF}
else
    lastal -l 3 -T 1 -e 20  -a 1 -b 1 -r 1 -q 1 ${BASE} ${FASTA2} >${MAF}
fi

maf-convert sam ${MAF}> ${SAM}


