#!/bin/bash
#ex. Forge/NewSkeletonWriter.sh ExampleRun.C NewName(.C) 

if [[ -z $SKFlat_WD ]]; then echo "source setup.sh needed."; exit 1; fi
if [[ -z $1 || -z $2 ]]; then echo "Missed either source or target."; exit 1; fi

TargetPath="$SKFlat_WD/Analyzers"
BaseNameSrc=$( echo $1 | rev | cut -d "/" -f1 | cut -d "." -f2 | rev )
BaseNameDest=$( echo $2 | rev | cut -d "/" -f1 | cut -d "." -f2 | rev )

if [[ ! -e ${TargetPath}/include/${BaseNameSrc}.h ]]; then echo "${BaseNameSrc}.h does not exist."; exit 1; fi
if [[ ! -e ${TargetPath}/src/${BaseNameSrc}.C ]]; then echo "${BaseNameSrc}.C does not exist."; exit 1; fi
if [[ -e ${TargetPath}/include/${BaseNameDest}.h ]]; then echo "${BaseNameDest}.h already exists."; exit 1; fi
if [[ -e ${TargetPath}/src/${BaseNameDest}.C ]]; then echo "${BaseNameDest}.C already exists."; exit 1; fi

cp ${TargetPath}/include/${BaseNameSrc}.h ${TargetPath}/include/${BaseNameDest}.h
cp ${TargetPath}/src/${BaseNameSrc}.C ${TargetPath}/src/${BaseNameDest}.C
sed -i "s/${BaseNameSrc}/${BaseNameDest}/g" ${TargetPath}/include/${BaseNameDest}.h 
sed -i "s/${BaseNameSrc}/${BaseNameDest}/g" ${TargetPath}/src/${BaseNameDest}.C 
sed -i "/my codes end/i #pragma link C++ class ${BaseNameDest}+;" ${TargetPath}/include/Analyzers_LinkDef.h 
