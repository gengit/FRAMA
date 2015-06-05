#!/usr/bin/env bash

# Arguments:
#   directory; filter; output

for i in $(find $1 -print | grep $2); do 
	echo $i; 
	dd if=$i of=$3 conv=notrunc oflag=append bs=5M
done
