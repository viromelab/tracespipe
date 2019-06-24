#!/bin/bash
#
index_f="../to_encrypt_data/index.txt";
echo "Please, enter the password to encrypt the files contained in ../to_encrypt_data directory: ";
read -s password
echo "$password" > key.txt
for file in ../to_encrypt_data/*
  do
  if [[ "$file" == "$index_f" ]];
  then
  echo "Jumping index file ($file)";
  else 
  echo "Encrypting $file ...";
  out_file=$(basename $file);
  cryfa -k key.txt $file > ../encrypted_data/$out_file.enc
  fi
  done
rm -f key.txt
echo "Done!";
#
