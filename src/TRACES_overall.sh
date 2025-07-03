#!/bin/bash
#
#USAGE: TRACES_overall.sh reftype pattern organ [maxdepth=1000]
#   RefType is one of viral, specific, cy, mtdna and defines the type of reference we are calculating stats for
#   Pattern defines the specific reference used
#       for cy and mtdna this should be cy and mt respectively
#       for specific it is the pattern used to find the reference in the DB
#       for virus it is the virus Name
#   Organ is the identifier library being processed
#   MaxDepth is the maximum depth to consider at any given position
#
#
#Define Global Variables
declare -A ValidRefTypes;
for type in cy mtdna specific viral; do
    ValidRefTypes[$type]=1
done 
declare -A RequiredPattern;
RequiredPattern["cy"]="cy"
RequiredPattern["mtdna"]="mt"
DEFAULT_MAX_DEPTH=1000;
#Pull the Ref Type, validate it, and set the directories containing relevant bed and bam, and output files
REFTYPE="$1";shift;
if [[ -z "${ValidRefTypes[$REFTYPE]}" ]]; then
    >&2 echo -e "[BUG] Attempt to call TRACES_overall.sh on an unrecognized type ($REFTYPE)\n" \
                "\t Valid types are: " "${!ValidRefTypes[@]}" \
                ;
    exit 1;
fi
BedDir="../output_data/TRACES_${REFTYPE}_bed"
BamDir="../output_data/TRACES_${REFTYPE}_alignments"
StatsDir="../output_data/TRACES_${REFTYPE}_statistics"
#Pull the pattern and validate it against the reftype
PATTERN="$1";shift;
if [[ -n "${RequiredPattern[$REFTYPE]}" && "$PATTERN" != "${RequiredPattern[$REFTYPE]}" ]]; then
    rq="${RequiredPattern[$REFTYPE]}"
    >&2 echo "[BUG WARNING] Call to TRACES_overall.sh $REFTYPE with pattern!='$rq'";
    PATTERN=$rq
fi
#Set relevant pref- and suffixes
BedPrefix="$PATTERN" #post-delimiter is constant across all ref types
AlnBamPrefix="$PATTERN" #post-delimiter is constant across all ref types
AlnBamSuffix="-$PATTERN" #pre-delimiter included as the suffix isn't present for cy and mtdna
StatsPrefix="$PATTERN" #post-delimiter is constant across all ref types
#Handle cases where the -fixes aren't pattern
case $REFTYPE in 
    viral)
        AlnBamPrefix="viral" 
        ;;
    specific)
        AlnBamPrefix="viral" 
        ;;
    cy)
        AlnBamSuffix=""
        ;;
    mtdna)
        AlnBamSuffix=""
        ;;
esac
#Pull the organ and then set the relevant bed, bam, and output files
ORGAN="$1"; shift
CovBedFile="$BedDir/$BedPrefix-coverage-$ORGAN.bed"
AlnBamFile="$BamDir/${AlnBamPrefix}_aligned_sorted-${ORGAN}${AlnBamSuffix}.bam"
DepthFile="$StatsDir/${StatsPrefix}-total-depth-coverage-${ORGAN}.txt"
ZeroBedFile="$BedDir/${BedPrefix}-zero-coverage-$ORGAN.bed"
BreadthFile="$StatsDir/${StatsPrefix}-total-horizontal-coverage-${ORGAN}.txt"
#Pull the Max Depth from the comand line options
MAX_DEPTH=$1; shift
if [ -z "$MAX_DEPTH" ]; then
    MAX_DEPTH=$DEFAULT_MAX_DEPTH;
fi
#
#Get the total size of the reference sequence, and the sequence ID for the first contig
#   Note: we assume all references are single contigs at this point
TOTAL_SIZE=$(tail -n 1 "$CovBedFile" | awk '{ print $3}');
acc=$(head -1 "$CovBedFile" | cut -f1);
bedEntry="$acc\t0\t$TOTAL_SIZE"
#Updated Depth calculation which ignores peaks, but more importantly ignores secondary/supplementary alignments
#Column 3 is the length of the genome (constructed in the bedEntry above
#Column 4 is the total bases covered
NORMALIZED_COVERAGE=$(samtools bedcov --max-depth "$MAX_DEPTH" <(echo -e "$bedEntry") "$AlnBamFile" | awk '{print $4/$3}')
printf "$PATTERN\t$ORGAN\t%0.4f\n" "$NORMALIZED_COVERAGE" > "$DepthFile"
#
ZERO_COVERAGE=$(awk '{sum += ($3-$2)} END {print sum}' "$ZeroBedFile");
if [ -z "$ZERO_COVERAGE" ]; then
    ZERO_COVERAGE=0;
fi
NORMALIZED_ZERO_COVERAGE=$(echo "scale=4; (100-(($ZERO_COVERAGE / $TOTAL_SIZE)*100))" | bc -l);
printf "$PATTERN\t$ORGAN\t%0.4f\n" "$NORMALIZED_ZERO_COVERAGE" > "$BreadthFile"
#
