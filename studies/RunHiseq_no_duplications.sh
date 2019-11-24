#!/bin/bash
#
echo "Running hiSeq ..."
cd hiseq_no_duplications/tracespipe/src/
./TRACESPipe.sh --run-meta --run-all-v-alig --remove-dup
cd ../../../
echo "Done!";
#
