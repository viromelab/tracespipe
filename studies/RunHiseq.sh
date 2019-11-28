#!/bin/bash
#
echo "Running hiSeq ..."
cd hiseq/tracespipe/src/
./TRACESPipe.sh --run-meta --run-all-v-alig 
cd ../../../
echo "Done!";
#
