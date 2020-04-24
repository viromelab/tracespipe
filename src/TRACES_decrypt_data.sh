#!/bin/bash
#
index_f="../encrypted_data/index.txt";
for file in ../encrypted_data/*
  do
  if [[ "$file" == "$index_f" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Jumping index file ($file) \e[0m";
  else 
  echo -e "\e[34m[TRACESPipe]\e[32m Decrypting $file \e[0m";
  out_file=$(basename $file);
  echo -e "\e[34m[TRACESPipe]\e[32m Please, enter the password: \e[0m";
  read -s password
  echo "$password" > key.txt
  cryfa -k key.txt -d $file > ../decrypted_data/$out_file.dec 2>> ../logs/Log-stderr-$ORGAN_T.txt;
  rm -f key.txt;
  fi
  done
echo -e "\e[34m[TRACESPipe]\e[32m Done! \e[0m";
#
