#!/bin/bash

########################################################################
#"runBkdMC" "runConv" "runFlip" "runFake" "runSigMC" "runData"
#declare -a List_runModes=("runBkdMC")
#declare -a List_runModes=("runBkdMC" "runData")
#declare -a List_runModes=("runSigMC")
#declare -a List_runModes=("runFake")
#declare -a List_runModes=("runFlip")
#declare -a List_runModes=("runBkdMC" "runFake" "runConv" "runFlip" "runData")
declare -a List_runModes=("runBkdMC" "runFake")
#declare -a List_runModes=("runBkdMC" "runFake" "runData")
#declare -a List_runModes=("runData" "runFake")
#declare -a List_runModes=("runBkdMC" "runData" "runFake")
#declare -a List_runModes=("runBkdMC" "runData" "runSigMC")

########################################################################
## RUN PARAMETERS

AnalysisCode="ControlPlots"
FinalState="OS2l" #"SS2l" "OS2l"
Skim="" #"SkimTree_SS2lOR3l"
Year="2017"
RunningMode="TopCR_OS2l,SystRun" #"SS2l,SB_SS2L,SystRun" 
#"SS2l,CFlip"
#"SS2l,TopBSrc" #"TopCR_OS2l" #"SS2l,SB_SS2L" #"TriLep,FkCR3l" #"TriLep,ConvCR" #"TriLep,ConvVar,ConvRun" #"FakeRun" #"ConvRun" #"FlipRun" #"SystRun"

NJobs=""
Memory=""
NJobMax="80"
NEvtMax="100000"
NSkipEvt=""
ReductionFactor="" #"10"
SpecificList=""
SpecificDir="Ver2" #"SB_SS2L/Ver3_2" #"SB_SS2L/Ver3" #"CFlip" #"Fk3lCR" #"FkCR3l" #"ConvCR" #"ConvVar"
runDebug="False"

# Options:
# List_runModes   - run on RD, BkdMC, SigMC, 
#                 - Opt(Mult): runData/runFake/runBkdMC/runSigMC
# FinalState      - LabelDir & Code options (flag)
#                 - Opt(Sigl): DoubleMuon... 
# Skim            - run on Skimmed data set.
# Year            - Opt(Sigl): 2016/2017/2018
# RunningMode     - userflag
# NJobs           - num. of jobs for each sample
#                 - Opt(Sigl): int - uniform job allocation / 
#                              ""  - job allocation from GetJobConfig.py
# Memory          - max. memory allocation per job.
# NJobMax         - max. num. of jobs permitted for user.
# NEvtMax         - Max. num of evts to run - available only for debug
#                 - Opt(Sigl): ""/int
# NSkipEvt        - Skip first designated num. of evts.
#                 - Opt(Sigl): ""/int
# ReductionFactor - run only Total#/ReductionFactor events
#                 - Opt(Sigl): ""/int
# runDebug        - Make a code for pilot running on master node
#                 - Opt(Sigl): True/False


########################################################################

if [[ ! -e CommandHist.txt ]]; then touch CommandHist.txt; echo "Command History">> CommandHist.txt; echo >> CommandHist.txt; fi;
date >> CommandHist.txt

for runMode in ${List_runModes[@]}; 
do

  ## OUTPUT PATH CONFIG
  CategoryLabel=""
  StateLabel=${FinalState}
  
  if [[ -z ${SKFlatV} ]]; then echo "source setup.sh needed. exit"; exit 1; fi
  if [[ ! -e "${SKFlat_WD}/Analyzers/include/${AnalysisCode}.h" ]]; then
    echo "Inappropriate AnalysisCode setting."; exit 1;
  fi
  if   [[ ${runMode} == "runData"  ]]; then CategoryLabel="Data";
  elif [[ ${runMode} == "runFake"  ]]; then CategoryLabel="Fake";
  elif [[ ${runMode} == "runConv"  ]]; then CategoryLabel="Conv";
  elif [[ ${runMode} == "runFlip"  ]]; then CategoryLabel="Flip";
  elif [[ ${runMode} == "runBkdMC" ]]; then CategoryLabel="BkdMC";
  elif [[ ${runMode} == "runSigMC" ]]; then CategoryLabel="SigMC";
  else echo "Inappropriate runData/runBkdMC/runFake/runSigMC settings."; exit 1; fi 

  OutputDir="${SKFlatOutputDir}/${SKFlatV}/${AnalysisCode}/${Year}/${StateLabel}/${SpecificDir}"
  if [[ ${runMode} == "runFake" || ${runMode} == "runConv" || ${runMode} == "runFlip" ]];
  then OutputDir+="/BkdMC/${CategoryLabel}"
  else OutputDir+="/${CategoryLabel}"; fi
  mkdir -pv ${OutputDir}
  
  
  #OPTION SETTING
  if [[ -z ${SpecificList} ]]; then
    List="${AnalysisCode}_${CategoryLabel}_${Year}";
    if [[ ${runDebug} == "True" ]]; then List="Debug_${CategoryLabel}_${Year}";
    elif [[ ${runMode} == "runData" || ${runMode} == "runFake" ]]; then List="${AnalysisCode}_Data_${Year}_${StateLabel}"; 
    elif [[ ${runMode} == "runFlip" ]]; then List="${AnalysisCode}_Flip_${Year}_${StateLabel}"; 
    else :; fi
  fi

  NSample=$( python GetJobConfig.py ${List} )
  #echo "python GetJobConfig.py ${List}"
  if [[ ${NSample} == "Error"* ]]; then echo "Sample List wrong for ${runMode}, skip."; continue; fi

  for (( it_s=0; it_s<${NSample}; it_s++ ))
  do
    Option=" -a ${AnalysisCode}"

    SampleName=$( python GetJobConfig.py ${List} ${it_s} | cut -d ' ' -f1 )  
    NJobs=$( python GetJobConfig.py ${List} ${it_s} | cut -d ' ' -f2 )
    DataPeriod="ALL"
    if [[ ${runMode} == "runData" || ${runMode} == "runFake" || ${runMode} == "runFlip" ]]; then DataPeriod=$( python GetJobConfig.py ${List} ${it_s} | cut -d ' ' -f3 ); fi

    if [[ ${NJobs} -gt 0 ]]; then :; 
    elif [[ ${NJobs} == "-" ]]; then NJobs=20; 
    else echo "NJob Config invalid for ${SampleName} in ${CategoryLabel} mode"; continue; fi
    if [[ ( ${runDebug} == "False" ) && ( -z ${ReductionFactor} ) ]]; then (( NJobs *= 2 )); fi

    TmpRunningMode=${RunningMode}; TmpSkim=${Skim};
    if [[ ${runMode} == "runFake" ]]; then TmpRunningMode+=",FakeRun"; fi
    if [[ ${runMode} == "runFlip" ]]; then TmpRunningMode+=",FlipRun"; TmpSkim=""; fi
    if [[ ${runMode} == "runConv" ]]; then TmpRunningMode+=",ConvRun"; fi

    Option+=" -i ${SampleName}" 
    Option+=" -n ${NJobs}" 
    Option+=" -y ${Year} -p ${DataPeriod}";
    if [[ ! -z ${TmpSkim} ]]; then Option+=" --skim ${TmpSkim}"; fi
    if [[ ! -z ${ReductionFactor} ]]; then Option+=" --reduction ${ReductionFactor}"; fi
    if [[ ! -z ${TmpRunningMode} ]]; then Option+=" --userflags ${TmpRunningMode}"; fi
    if [[ ! -z ${Memory} ]]; then Option+=" --memory ${Memory}"; fi
    if [[ ! -z ${NJobMax} ]]; then Option+=" --nmax ${NJobMax}"; fi
    if [[ ${runDebug} == "True" ]]; then Option+=" --no_exec"; fi
    Option+=" -o ${OutputDir}"
    SKFlat.py ${Option} &
    echo "SKFlat.py ${Option}" >> CommandHist.txt

    if [[ ${runDebug} == "True" ]]; then 
      sleep 4s
      DirName=$( ls -rt ${SKFlatRunlogDir} | tail -1 ) 
      if [[ ${DirName} == "www"* ]]; then DirName= $( ls -rt ${SKFlatRunlogDir} | tail -2 | head -1 ); fi
      PeriodLabel=""
      if [[ ${runMode} == "runData" || ${runMode} == "runFake" || ${runMode} == "runFlip" ]]; then PeriodLabel="_period${DataPeriod}"; fi
      OrigPath=${SKFlatRunlogDir}/${DirName}/${SampleName}${PeriodLabel}
      DebugPath=${SKFlat_WD}/Debug_${AnalysisCode}_${SampleName}${PeriodLabel}
      mkdir -p ${DebugPath}
      mv ${OrigPath}/run_0.C ${DebugPath}/run_0.C
      sed -i "1 i R__LOAD_LIBRARY(${SKFlat_WD}/lib/libDataFormats.so)" ${DebugPath}/run_0.C
      sed -i "1 i R__LOAD_LIBRARY(${SKFlat_WD}/lib/libAnalyzers.so)" ${DebugPath}/run_0.C
      if [[ ${NEvtMax} -gt 0 ]]; then sed -i "/m.LogEvery/a \ \ m.MaxEvent\ =\ ${NEvtMax};" ${DebugPath}/run_0.C; fi
      if [[ ${NSkipEvt} -gt 0 ]]; then sed -i "/m.LogEvery/a \ \ m.NSkipEvent\ =\ ${NSkipEvt};" ${DebugPath}/run_0.C; fi
    fi
  done
done

echo >> CommandHist.txt
