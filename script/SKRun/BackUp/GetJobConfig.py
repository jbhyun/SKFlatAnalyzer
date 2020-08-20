import os, sys


List_NewAnalyzer_SigMC_2017 = [
  ["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1], 
  ["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA85", 1],
  ["TestUnFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1],
  ["TestUnFiltered_TTToHcToWA_AToMuMu_MHc160_MA85", 1] 
]

List_NewAnalyzer_Data_2017_DoubleMuon = [
  ["DoubleMuon", "-", "B"],
  ["DoubleMuon", "-", "C"],
  ["DoubleMuon", "-", "D"],
  ["DoubleMuon", "-", "E"],
  ["DoubleMuon", "-", "F"],
]


List_Debug_SigMC_2017 = [
  ["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1], 
]

List_Debug_BkdMC_2017 = [
  ["ttZToLLNuNu", 1], 
]

List_Debug_Data_2017 = [
  ["DoubleMuon", 1, "B"], 
]


Dict_JobConfig = {
 "NewAnalyzer_SigMC_2017": List_NewAnalyzer_SigMC_2017,
 "Debug_SigMC_2017": List_Debug_SigMC_2017,
 "Debug_BkdMC_2017": List_Debug_BkdMC_2017,
 "Debug_Data_2017": List_Debug_Data_2017,
}


List = Dict_JobConfig.get(sys.argv[1])
if not List is None:
  if len(sys.argv) < 3:
    print(len(List))
  else:
    if "MC" in sys.argv[1]:
      print(str(List[int(sys.argv[2])][0])+" "+str(List[int(sys.argv[2])][1]))
    elif "Data" in sys.argv[1] or "Fake" in sys.argv[1]:
      print(str(List[int(sys.argv[2])][0])+" "+str(List[int(sys.argv[2])][1])+" "+str(List[int(sys.argv[2])][2]))
    else:
      print("invalid runMode option")
else:
  print("Error : No such list")
