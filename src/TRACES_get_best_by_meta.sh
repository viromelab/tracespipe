#!/bin/bash
ORGAN=$1; shift
Virus=$1; shift
MetaFile=$1; shift
#
awk -v Virus="$Virus" -F '\\t' '
    BEGIN{OFS="\t";bestSim="-";bestAcc="-"}
    (ARGIND == 1){if($2==Virus){InSet[$1]=1};next}
    {
        sub("NC_","NC-",$4);
        split($4,a,"_");
        acc = a[1];
        sub("NC-","NC_",acc);
        sim=$3;
    }
    (InSet[acc] && sim > bestSim+0){ bestSim=sim; bestAcc=acc;}
    END {print bestSim,bestAcc}
' "$MetaFile" "top-$ORGAN.csv"
