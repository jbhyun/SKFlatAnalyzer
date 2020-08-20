#!/bin/bash

if [[ -z $SKFlatRunlogDir ]]; then echo "source setup.sh needed."; exit 1; fi;
rm -r $SKFlatRunlogDir/*
