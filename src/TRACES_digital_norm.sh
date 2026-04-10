#!/bin/bash

function ConstructOutPath {
    path=$1; shift
    echo "$(dirname "$path")/DN-$(basename "$path")"
}

function CreateLink {
    path=$1; shift
    ln -s "$(readlink -f "$path")" "$(ConstructOutPath "$path")";
}

#Taken from stackoverflow question 1527049
# Answer from user Nicholas Sushkin
function JoinBy {
	local d=${1-} f=${2-}
	if shift 2; then
		printf %s "$f" "${@/#/$d}"
	fi
}

function InterleaveFastQ {
    fq1=$1; shift
    fq2=$1; shift
    awk -v fstr="$fq1,$fq2" '
        BEGIN{
            n=split(fstr,files,",");
            bCont=1;
            while(bCont){
                for(fi=1;fi<=n && bCont;fi++){
                    for(i=1;i<=4 && bCont;i++){
                        if(!(getline < (files[fi]))){
                            if(i > 1 || fi > 1){
                                print "[WARNING] Incomplete Fastq entry or unbalanced input when interleaving fastqs" > "/dev/stderr"
                            }
                            bCont = 0;
                            continue
                        };
                        print
                    }
                }
            }
        }
    '
}

#BBMap assumes paired reads are the same length in bytes!
#As preprocessing can violate that assumption, this leads to different numbers
# of reads in the same buffer size, leading to failed assertions
#As a workaround, we can interleave the input. In the default 2 pass mode,
# bbnorm creates temporary paired files, and reads them in the same way leading
# to the same assertion errors. One could run one pass mode, re-interleave the
# input and run a second pass but I'm not certain this is the same as 1 2-pass
# run, so instead we'll do just a single pass
function InterleavedBBNorm {
    in1=$1; shift
    in2=$1; shift
    extraArg=$1; shift
    threads=$1; shift
    out1=$1; shift
    out2=$1; shift
    tmpFile="tmp.$(basename "$in1" .fq)-$(basename "$in2" .fq)-interleaved.fq"
    rm -f "$tmpFile;"
    if ! InterleaveFastQ "$in1" "$in2" > "$tmpFile"; then
        >&2 echo "[ERROR] Interleaved BBNorm - Failure to interleave input fqs"
        rm -f "$tmpFile"
        return 1;
    fi
    if ! bbnorm.sh -in="$tmpFile" interleaved=true passes=1 "$extraArg" \
        threads="$maxThreads" out="$(ConstructOutPath "$fwdPath")" \
        out2="$(ConstructOutPath "$revPath")";
    then
        >&2 echo "[ERROR] Interleaved BBNorm - BBNorm Failure - Check $tmpFile"
        return 1;
    fi
    #Cleanup interleaved input
    rm -f "$tmpFile"
}

if [ "$#" -lt 6 ]; then
    >&2 echo "Usage: $(basename "$0") maxMem memFactor threads PathTo/rawFwd.gz fwd.fq rev.fq [ unpaired1.fq ... ]";
    >&2 echo -e "\tFor the fwd and rev paths creates a file in the same directory with the prefix 'DN-'\n" \
                "\t\tThat file will be a link if no normalization occurs\n" \
                "\tunpaired reads are used for kmer hashing, but are not included in the output\n" \
                "\t\tNo new files or links will be generated for them\n";
    exit 1;
fi

maxMem=$1;shift
memFactor=$1;shift
maxThreads=$1;shift
rawFwdGZPath=$1;shift
fwdPath=$1;shift
revPath=$1;shift
extraArg="";
if [ "$#" -gt 0 ]; then
    extraArg="extra=$(JoinBy "," "$@")"
fi

bFilter=0;
#Check if digital normalization is required
if [[ "$maxMem" -gt 0 && "$memFactor" -gt 0 ]]; then
    bFilter=$(stat -Lc "%s*$memFactor > $maxMem*1024^3" "$rawFwdGZPath" | bc);
fi

if [ "$bFilter" -eq 1 ]; then
    ##Ideal way if BBNorm gets fixed
    #if ! bbnorm.sh in="$fwdPath" in2="$revPath" "$extraArg" \
    #    threads="$maxThreads" out="$(ConstructOutPath "$fwdPath")" \
    #    out2="$(ConstructOutPath "$revPath")"
    ##Attempt digital normalization
    if ! InterleavedBBNorm "$fwdPath" "$revPath" "$extraArg" "$threads"; then
        bFilter=0;
        >&2 echo -e "\e[33mWARNING\e[0m: Digital Normalization Failed";
    fi;
fi

if [ "$bFilter" -eq 0 ]; then
    CreateLink "$fwdPath"
    CreateLink "$revPath"
fi
