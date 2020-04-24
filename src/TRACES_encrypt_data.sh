#!/bin/bash
#
index_f="../to_encrypt_data/index.txt";
echo -e "\e[34m[TRACESPipe]\e[32m Please, enter the password to encrypt the files contained in ../to_encrypt_data directory: \e[0m";
read -s password
echo "$password" > key.txt
for file in ../to_encrypt_data/*
  do
  if [[ "$file" == "$index_f" ]];
  then
  echo -e "\e[34m[TRACESPipe]\e[32m Jumping index file ($file) \e[0m";
  else 
  echo -e "\e[34m[TRACESPipe]\e[32m Encrypting $file \e[0m";
  out_file=$(basename $file);
  cryfa -k key.txt $file > ../encrypted_data/$out_file.enc 2>> ../logs/Log-stderr-$ORGAN_T.txt;
  fi
  done
rm -f key.txt
echo -e "\e[34m[TRACESPipe]\e[32m Done! \e[0m";
#
