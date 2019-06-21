#!/bin/bash
#
MAX_THREADS=`cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l`;
#
echo -e "\e[32mMaximum number of threads of this machine: $MAX_THREADS \e[0m";
#
