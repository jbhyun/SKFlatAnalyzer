#!/bin/bash

List=$1

if [[ ! -e $List ]]; then echo "No List"; exit 1; fi

while read line
do
  kill -9 $line
done<$List
