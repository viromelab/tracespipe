#!/bin/bash
ORGAN=$1; shift
Virus=$1; shift
MetaFile=$1; shift
#
awk -v Virus="$Virus" -F '\\t' '
    BEGIN{OFS="\t";bestSim="-";bestAcc="-";bestLen="-"}
    (ARGIND == 1){if($2==Virus){InSet[$1]=1};next}
    {
        sub("NC_","NC-",$4);
        split($4,a,"_");
        acc = a[1];
        sub("NC-","NC_",acc);
        sim=$3;
        len=$2;
    }
    (InSet[acc] && sim > bestSim+0){ bestSim=sim; bestAcc=acc; bestLen=len;}
    END {print bestSim,bestAcc,bestLen}
' "$MetaFile" "top-$ORGAN.csv"
