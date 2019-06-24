#!/bin/bash
#
index_f="../encrypted_data/index.txt";
for file in ../encrypted_data/*
  do
  if [[ "$file" == "$index_f" ]];
  then
  echo "Jumping index file ($file)";
  else 
  echo "Decrypting $file ...";
  out_file=$(basename $file);
  echo "Please, enter the password: ";
  read -s password
  echo "$password" > key.txt
  cryfa -k key.txt -d $file > ../decrypted_data/$out_file.dec
  rm -f key.txt;
  fi
  done
echo "Done!";
#
