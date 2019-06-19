#!/bin/bash
#
for file in V*fastq.gz
  do
  echo "Running $file ...";
  cryfa -v -k pass.txt $file > $file-cryfed
  echo "Done!";
  done
#
