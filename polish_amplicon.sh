#!/bin/bash

#  Wed Feb 24 20:17:11 2016
#  kpalin@merit.ltdk.helsinki.fi

# For debugging
#set -o verbose 

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail



THREADS=16

IDENTITY=0.6


usage()  {
    echo -e "usage:\n$0 [-c WORKDIR] input.fasta

-c WORKDIR  Continue work on this directory
-h          Show this message and exit." >&2
    exit 1
    
}

CWD=$(dirname $(readlink -f $0)) 



while getopts  "hc:" flag
do
    case "$flag" in 
	c)
	TEMPDIR="$OPTARG"
	;;
	h|*)
	usage
	;;
    esac
done
shift $((OPTIND-1)); OPTIND=1


INFASTA=$(readlink -f $1)
echo "Using input:" $INFASTA


if [ ! -e ${TEMPDIR:-XXX} ];
then

    # Make temporary directory which will be cleaned at script exit
    TEMPDIR=$(mktemp --directory --tmpdir=$PWD AMPCORRECT.$(date -Is).XXXX)
#function _cleanup {
#    echo $TEMPDIR
#}
#trap _cleanup EXIT
fi

cd $TEMPDIR

test -s basecalls_sorted.fasta || sizeseq -sequences "${INFASTA}" -descending -outseq basecalls_sorted.fasta
test -s clusters.fasta || sumaclust -e -p ${THREADS} -t ${IDENTITY} -s None  basecalls_sorted.fasta  >clusters.fasta
test -s fasta_clusters.lst || python ${CWD}/split_fasta_by_cluster.py clusters.fasta > fasta_clusters.lst


# Only correct clusters with at least three members
if [ ! -s poa.batch ];
then
    for i in $(awk -v ORS=" " '$2>=3 {print $1;}' fasta_clusters.lst )
    do
	echo -read_fasta $i -clustal $(basename $i .fasta).clustal  $CWD/poa-blosum80.mat
    done|slargs -N -J poaV2 poa -tolower -do_global -hb >poa.batch
fi

SJOBstr="Dummy"

SJOBstr=$(sbatch --parsable <poa.batch)

echo $SJOBstr

cat >summarise.sh <<EOF 
#!/bin/bash
#SBATCH -J corrected
#SBATCH --dependency afterany:${SJOBstr}
#SBATCH --cpus-per-task=4

#  Wed Feb 24 20:17:11 2016
#  kpalin@merit.ltdk.helsinki.fi

# For debugging
#set -o verbose 

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail


awk  '\$2>=3 {sub(/[.]fasta\$/,".clustal",\$1);print \$1}' fasta_clusters.lst |xargs python $CWD/clustal2fasta.py  >correctedAmplicons.fasta
awk  '\$2<3 {print \$1}' fasta_clusters.lst |xargs cat >>correctedAmplicons.fasta

bwa mem -t 4  /mnt/cg8/reference-genomes/hs37d5_viral_bwa0.7.12/hs37d5_viral.fa correctedAmplicons.fasta | samtools view -Sb - | samtools sort - correctedAmplicons.sorted 

samtools index correctedAmplicons.sorted.bam

samtools stats correctedAmplicons.sorted.bam|grep -E '^(#|SN)'

EOF

sbatch <summarise.sh
