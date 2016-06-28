#!/usr/bin/env bash
if [ "$1" == "${1#/}" ] ; then inf=$(pwd)/$1 ; else inf=$1 ; fi
tmp=/tmp/polyar$$
ssh gen101 "( cd /home/lamarck/mbens/app/src/polyar/ ; java -classpath ./ polyar -i $inf -o $tmp ; cat $tmp )"
rm $tmp
