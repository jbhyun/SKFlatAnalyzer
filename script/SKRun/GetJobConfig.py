import os, sys


List_NewAnalyzer_SigMC_2017 = [
  ["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1], ["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA85", 1], ["TestUnFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1], ["TestUnFiltered_TTToHcToWA_AToMuMu_MHc160_MA85", 1] 
]

List_DiLepValid_Data_2017_DiElectron = [
  ["SingleElectron", 9, "B"], ["SingleElectron", 21, "C"], ["SingleElectron", 8, "D"], ["SingleElectron", 16, "E"], ["SingleElectron", 21, "F"],
]
List_DiLepValid_Data_2017_DiMuon = [
  ["SingleMuon", 13, "B"], ["SingleMuon", 15, "C"], ["SingleMuon", 7, "D"], ["SingleMuon", 15, "E"], ["SingleMuon", 24, "F"],
]
List_DiLepValid_Data_2018_DiElectron = [
  ["EGamma", 18, "A"], ["EGamma", 8, "B"], ["EGamma", 8, "C"], ["EGamma", 41, "D"],
]
List_DiLepValid_Data_2018_DiMuon = [
  ["SingleMuon", 18, "A"], ["SingleMuon", 9, "B"], ["SingleMuon", 9, "C"], ["SingleMuon", 39, "D"],
]
List_DiLepValid_BkdMC_2017 = [
  ["DYJets", 21], ["TTLL_powheg", 13], ["WZ_pythia", 1], ["WW_pythia", 2], ["ZZ_pythia", 1], ["SingleTop_tW_antitop_NoFullyHad", 2], ["SingleTop_tW_top_NoFullyHad", 2],
]

List_BTagEff_MC_BkdMC_2017 = [["TTLJ_powheg", 8], ["TTLL_powheg", 13], ["DYJets", 21]]
List_BTagEff_MC_BkdMC_2018 = [["TTLJ_powheg", 50], ["TTLL_powheg", 13]]

List_TestRun_Data_2016_TriMu = [
  ["DoubleMuon", 15, "B_ver2"], ["DoubleMuon", 7, "C"], ["DoubleMuon", 11, "D"], ["DoubleMuon", 11, "E"], ["DoubleMuon", 8, "F"], ["DoubleMuon", 16, "G"], ["DoubleMuon", 14, "H"],
]
List_TestRun_Data_2016_ElDiMu = [
  ["MuonEG", 15, "B_ver2"], ["MuonEG", 7, "C"], ["MuonEG", 11, "D"], ["MuonEG", 11, "E"], ["MuonEG", 8, "F"], ["MuonEG", 16, "G"], ["MuonEG", 14, "H"],
]
List_TestRun_Data_2017_TriMu = [
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
]
List_TestRun_Data_2017_ElDiMu = [
  ["MuonEG", 13, "B"], ["MuonEG", 15, "C"], ["MuonEG", 7, "D"], ["MuonEG", 15, "E"], ["MuonEG", 24, "F"],
]
List_TestRun_Data_2018_TriMu = [
  ["DoubleMuon", 18, "A"], ["DoubleMuon", 9, "B"], ["DoubleMuon", 9, "C"], ["DoubleMuon", 39, "D"],
]
List_TestRun_Data_2018_ElDiMu = [
  ["MuonEG", 18, "A"], ["MuonEG", 9, "B"], ["MuonEG", 9, "C"], ["MuonEG", 39, "D"],
]


List_HNTopFeas_BkdMC_2016 = [
  ["DYJets", 25], ["TTLL_powheg", 13],
  ["ZZTo4L_powheg", 20], ["ggHToZZTo4L", 2], ["VBF_HToZZTo4L", 1], ["ttHToNonbb", 2], ["ttZToLLNuNu", 4], ["WWZ", 1], ["WZZ", 1], ["ZZZ", 1],
  ["WZTo3LNu_powheg", 5], ["ttWToLNu", 2], ["WWW", 1], #Only Trilep  
  ["WJets_MG", 13], ["TTLJ_powheg", 25], #Only SS2l
  ["WGToLNuG", 5], ["ZGTo2LG", 5], ["TTG", 5], ["TG", 1], ["VHToNonbb", 1], ["WpWp_EWK", 1], ["WpWp_QCD", 1],
]
List_HNTopFeas_BkdMC_2017 = [
  ["DYJets", 25], ["TTLL_powheg", 13],
  ["ZZTo4L_powheg", 3], ["ggHToZZTo4L", 2], ["VBF_HToZZTo4L", 1], ["ttHToNonbb", 2], ["ttZToLLNuNu", 4], ["WWZ", 1], ["WZZ", 1], ["ZZZ", 1],
  ["WZTo3LNu_powheg", 1], ["ttWToLNu", 2], ["WWW", 1], #Only Trilep  
  ["WJets_MG", 7], ["TTLJ_powheg", 7], #Only SS2l
  ["WGToLNuG_01J", 5], ["ZGToLLG_01J", 5], ["TTG", 2], ["TG", 1], ["VHToNonbb", 1], ["WpWp_EWK", 1], ["WpWp_QCD", 1],
]
List_HNTopFeas_BkdMC_2018 = [
  ["DYJets", 25], ["TTLL_powheg", 13],
  ["ZZTo4L_powheg", 20], ["ggHToZZTo4L", 2], ["VBF_HToZZTo4L", 1], ["ttHToNonbb", 2], ["ttZToLLNuNu", 4], ["WWZ", 1], ["WZZ", 1], ["ZZZ", 1],
  ["WZTo3LNu_powheg", 1], ["ttWToLNu", 2], ["WWW", 1], #Only Trilep  
  ["WJets_MG", 13], ["TTLJ_powheg", 25], #Only SS2l
  ["WGToLNuG_01J", 5], ["ZGToLLG_01J", 5], ["TTG", 2], ["TG", 1], ["VHToNonbb", 1], ["WpWp_EWK", 1], ["WpWp_QCD", 1],
]


List_SkimTree_SS2lOR3l_Data_2016 = [
  ["DoubleMuon", 15, "B_ver2"], ["DoubleMuon", 7, "C"], ["DoubleMuon", 11, "D"], ["DoubleMuon", 11, "E"], ["DoubleMuon", 8, "F"], ["DoubleMuon", 16, "G"], ["DoubleMuon", 14, "H"],
  ["MuonEG", 15, "B_ver2"], ["MuonEG", 7, "C"], ["MuonEG", 11, "D"], ["MuonEG", 11, "E"], ["MuonEG", 8, "F"], ["MuonEG", 16, "G"], ["MuonEG", 14, "H"],
  ["DoubleEG", 15, "B_ver2"], ["DoubleEG", 7, "C"], ["DoubleEG", 11, "D"], ["DoubleEG", 11, "E"], ["DoubleEG", 8, "F"], ["DoubleEG", 16, "G"], ["DoubleEG", 14, "H"],
]
List_SkimTree_SS2lOR3l_Data_2017 = [
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
  ["MuonEG", 13, "B"], ["MuonEG", 15, "C"], ["MuonEG", 7, "D"], ["MuonEG", 15, "E"], ["MuonEG", 24, "F"],
  ["DoubleEG", 13, "B"], ["DoubleEG", 15, "C"], ["DoubleEG", 7, "D"], ["DoubleEG", 15, "E"], ["DoubleEG", 24, "F"],
]
List_SkimTree_SS2lOR3l_Data_2018 = [
  ["DoubleMuon", 18, "A"], ["DoubleMuon", 9, "B"], ["DoubleMuon", 9, "C"], ["DoubleMuon", 39, "D"],
  ["MuonEG", 18, "A"], ["MuonEG", 9, "B"], ["MuonEG", 9, "C"], ["MuonEG", 39, "D"],
  ["EGamma", 18, "A"], ["EGamma", 9, "B"], ["EGamma", 9, "C"], ["EGamma", 39, "D"],
]


List_SkimRateCheck_BkdMC_2016 = [
  ["DYJets", 25], ["TTLL_powheg", 13], ["WJets_MG", 13], ["TTLJ_powheg", 25], #Only SS2l
]

List_TrigCheck_BkdMC_2016 = [
  ["WZTo3LNu_powheg", 5], ["ZZTo4L_powheg", 20],
]

List_TrigCheck_BkdMC_2017 = [
  ["TTLL_powheg", 13], #Only SS2l
]


List_MCPUDist17_BkdMC_2017 = [
  ["TT_TTobNMu_SS2L_LO_MN20", 1],
  ["TT_TTobNMu_SS2L_LO_MN50", 1],
  ["TT_TTobNMu_SS2L_LO_MN100", 1],
]

List_GenMatchingValid_BkdMC_2016 = [
  ["ZGTo2LG", 5], ["TTLJ_powheg", 25],
]

List_GenMatchingValid_BkdMC_2017 = [
  ["ZGToLLG_01J", 5], ["TTLJ_powheg", 7],
]

List_GenMatchingValid_BkdMC_2018 = [
  ["ZGToLLG_01J", 5], ["TTLJ_powheg", 25],
]

List_GenMatchingValid_SigMC_2017 = [
  ["TT_TTobNMu_SS2L_LO_MN20", 1], ["TT_TTobNMu_SS2L_LO_MN50", 1], ["TT_TTobNMu_SS2L_LO_MN100", 1],
]






List_Debug_SigMC_2016 = [["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1]]
List_Debug_BkdMC_2016 = [
  ["TTLL_powheg", 1]
  #["ZGTo2LG", 1]
]
List_Debug_Data_2016 = [["MuonEG", 1, "B_ver2"]]
List_Debug_SigMC_2017 = [["TT_TTobNMu_SS2L_LO_MN50", 1]]
List_Debug_BkdMC_2017 = [["ttZToLLNuNu", 1]]
List_Debug_Data_2017 = [["SingleMuon", 1, "B"]]
List_Debug_SigMC_2018 = [["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1]]
List_Debug_BkdMC_2018 = [["ZZTo4L_powheg", 1]]
List_Debug_Data_2018 = [["SingleMuon", 1, "B"]]



Dict_JobConfig = {
 "GenMatchingValid_BkdMC_2016": List_GenMatchingValid_BkdMC_2016,
 "GenMatchingValid_BkdMC_2017": List_GenMatchingValid_BkdMC_2017,
 "GenMatchingValid_BkdMC_2018": List_GenMatchingValid_BkdMC_2018,
 "GenMatchingValid_SigMC_2017": List_GenMatchingValid_SigMC_2017,
 "TrigCheck_Data_2017_DiMuon": List_TestRun_Data_2017_TriMu,
 "TrigCheck_BkdMC_2016": List_TrigCheck_BkdMC_2016,
 "TrigCheck_BkdMC_2017": List_TrigCheck_BkdMC_2017,
 "SkimRateCheck_BkdMC_2016": List_SkimRateCheck_BkdMC_2016,
 "SkimRateCheck_Data_2016_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2016,
 "SkimTree_SS2lOR3l_BkdMC_2016": List_HNTopFeas_BkdMC_2016,
 "SkimTree_SS2lOR3l_BkdMC_2017": List_HNTopFeas_BkdMC_2017,
 "SkimTree_SS2lOR3l_BkdMC_2018": List_HNTopFeas_BkdMC_2018,
 "SkimTree_SS2lOR3l_Data_2016_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2016,
 "SkimTree_SS2lOR3l_Data_2017_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2017,
 "SkimTree_SS2lOR3l_Data_2018_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2018,
 "SyncYield_Data_2016_ElDiMu": List_TestRun_Data_2016_ElDiMu,
 "SyncYield_Data_2016_TriMu": List_TestRun_Data_2016_TriMu,
 "MCPUDist17_BkdMC_2017" : List_MCPUDist17_BkdMC_2017,
 "HNTopFeas_BkdMC_2016": List_HNTopFeas_BkdMC_2016,
 "HNTopFeas_BkdMC_2017": List_HNTopFeas_BkdMC_2017,
 "HNTopFeas_BkdMC_2018": List_HNTopFeas_BkdMC_2018,
 "TestRun_Data_2016_TriMu": List_TestRun_Data_2016_TriMu,
 "TestRun_Data_2017_TriMu": List_TestRun_Data_2017_TriMu,
 "TestRun_Data_2018_TriMu": List_TestRun_Data_2018_TriMu,
 "TestRun_Data_2016_ElDiMu": List_TestRun_Data_2016_ElDiMu,
 "TestRun_Data_2017_ElDiMu": List_TestRun_Data_2017_ElDiMu,
 "TestRun_Data_2018_ElDiMu": List_TestRun_Data_2018_ElDiMu,
 "DiLepValid_Data_2017_DiElectron": List_DiLepValid_Data_2017_DiElectron,
 "DiLepValid_Data_2017_DiMuon": List_DiLepValid_Data_2017_DiMuon,
 "DiLepValid_Data_2017_ElectronMuon": List_DiLepValid_Data_2017_DiMuon,
 "DiLepValid_Data_2018_DiElectron": List_DiLepValid_Data_2018_DiElectron,
 "DiLepValid_Data_2018_DiMuon": List_DiLepValid_Data_2018_DiMuon,
 "DiLepValid_Data_2018_ElectronMuon": List_DiLepValid_Data_2018_DiMuon,
 "DiLepValid_BkdMC_2016": List_DiLepValid_BkdMC_2017,
 "DiLepValid_BkdMC_2017": List_DiLepValid_BkdMC_2017,
 "DiLepValid_BkdMC_2018": List_DiLepValid_BkdMC_2017,
 "BTagEff_MC_BkdMC_2017" : List_BTagEff_MC_BkdMC_2017,
 "BTagEff_MC_BkdMC_2018" : List_BTagEff_MC_BkdMC_2018,
 "NewAnalyzer_SigMC_2017": List_NewAnalyzer_SigMC_2017,
 "Debug_SigMC_2016": List_Debug_SigMC_2016,
 "Debug_BkdMC_2016": List_Debug_BkdMC_2016,
 "Debug_Data_2016": List_Debug_Data_2016,
 "Debug_SigMC_2017": List_Debug_SigMC_2017,
 "Debug_BkdMC_2017": List_Debug_BkdMC_2017,
 "Debug_Data_2017": List_Debug_Data_2017,
 "Debug_SigMC_2018": List_Debug_SigMC_2018,
 "Debug_BkdMC_2018": List_Debug_BkdMC_2018,
 "Debug_Data_2018": List_Debug_Data_2018,
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
