import os, sys

List_MC_AllBkd3l_2017 = [
  ["DYJets", 25], ["TTLL_powheg", 13],
  ["WZTo3LNu_powheg", 1], ["ZZTo4L_powheg", 3], ["ttWToLNu", 2], ["ttZToLLNuNu", 4], ["ttHToNonbb", 2],
  ["ggHToZZTo4L", 2], ["VBF_HToZZTo4L", 1], ["VHToNonbb", 1], ["WWZ", 1], ["WZZ", 1], ["ZZZ", 1],["WWW", 1], 
  ["ZGToLLG_01J", 5], ["TTG", 2],
]

List_MC_CF2l_2017 = [
  #["DYJets", 40],
  ["DYJets", 25], ["TTLL_powheg", 13], ["SingleTop_tW_antitop_NoFullyHad", 2], ["SingleTop_tW_top_NoFullyHad", 2], 
]

List_MC_CFMeas_2017 = [
# ["DYJets_MG_HT-40to70", 10], ["DYJets_MG_HT-70to100", 10], ["DYJets_MG_HT-100to200", 10], ["DYJets_MG_HT-200to400", 10], ["DYJets_MG_HT-400to600", 10], ["DYJets_MG_HT-600to800", 10], ["DYJets_MG_HT-800to1200", 5], ["DYJets_MG_HT-1200to2500", 5], ["DYJets_MG_HT-2500toInf", 5], 
 ["DYJets", 25],
]


List_MC_Cv2l_2017 = [
  #["ZGToLLG_01J", 5],
  ["WJets_MG", 13], ["DYJets", 25], ["TTLJ_powheg", 25], ["TTLL_powheg", 13], 
  ["ZGToLLG_01J", 5], ["WGToLNuG_01J", 5], ["TTG", 2], ["TG", 1], 
  ["SingleTop_tW_antitop_NoFullyHad", 2], ["SingleTop_tW_top_NoFullyHad", 2], ["SingleTop_sch_Lep", 5], ["SingleTop_tch_antitop_Incl", 5], ["SingleTop_tch_top_Incl", 5],
]

List_MC_Pr2l_2017 = [
  ["WZTo3LNu_powheg", 1], ["ZZTo4L_powheg", 3], ["ttWToLNu", 2], ["ttHToNonbb", 2], ["ttZToLLNuNu", 4],
  ["ggHToZZTo4L", 2], ["VBF_HToZZTo4L", 1], ["VHToNonbb", 1], ["WWW", 1], ["WWZ", 1], ["WZZ", 1], ["ZZZ", 1], ["WpWp_EWK", 1], ["WpWp_QCD", 1],
]

List_MC_Cv3l_2017 = [
  ["ZGToLLG_01J", 5], ["TTG", 2], ["DYJets", 25], 
]

List_MC_Fk3l_2017 = [
  #["DYJets", 40],
  ["DYJets", 25], ["TTLL_powheg", 13],
]

List_MC_NonFk3l_2017 = [
  ["WZTo3LNu_powheg", 1], ["ZZTo4L_powheg", 3], ["ttWToLNu", 5], ["ttZToLLNuNu", 5], ["ttHToNonbb", 5],
  ["ggHToZZTo4L", 2], ["VBF_HToZZTo4L", 1], ["VHToNonbb", 1], ["WWZ", 1], ["WZZ", 1], ["ZZZ", 1],["WWW", 1], 
#  ["ZGToLLG_01J", 5], ["TTG", 2],
]

List_Data_2LPD_2017 = [
  ["DoubleEG", 9, "B"], ["DoubleEG", 21, "C"], ["DoubleEG", 8, "D"], ["DoubleEG", 16, "E"], ["DoubleEG", 21, "F"],
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
  ["MuonEG", 13, "B"], ["MuonEG", 15, "C"], ["MuonEG", 7, "D"], ["MuonEG", 15, "E"], ["MuonEG", 24, "F"],
]

List_Data_SF2L_2017 = [
  ["DoubleEG", 9, "B"], ["DoubleEG", 21, "C"], ["DoubleEG", 8, "D"], ["DoubleEG", 16, "E"], ["DoubleEG", 21, "F"],
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
]

List_Data_DiEl_2017 = [
  ["DoubleEG", 9, "B"], ["DoubleEG", 21, "C"], ["DoubleEG", 8, "D"], ["DoubleEG", 16, "E"], ["DoubleEG", 21, "F"],
]


########################################################################################
List_NewAnalyzer_SigMC_2017 = [
  ["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1], ["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA85", 1], ["TestUnFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1], ["TestUnFiltered_TTToHcToWA_AToMuMu_MHc160_MA85", 1] 
]

List_FakeRateMeas_Data_2017_Mu = [
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
]
List_FakeRateMeas_Data_2017_El = [
  ["SingleElectron", 9, "B"], ["SingleElectron", 21, "C"], ["SingleElectron", 8, "D"], ["SingleElectron", 16, "E"], ["SingleElectron", 21, "F"],
]
List_FakeRateMeas_BkdMC_2017 = [
  ["WJets_MG", 13],  ["DYJets", 21], ["DYJets10to50_MG", 13],
  ["TTLJ_powheg", 25],  ["TTLL_powheg", 13], ["WZ_pythia", 1], ["WW_pythia", 2], ["ZZ_pythia", 1], ["SingleTop_tW_antitop_NoFullyHad", 2], ["SingleTop_tW_top_NoFullyHad", 2],
  ["SingleTop_sch_Lep", 5], ["SingleTop_tch_antitop_Incl", 5], ["SingleTop_tch_top_Incl", 5],
#  ["QCD_Pt-1000toInf_MuEnrichedPt5", 5], ["QCD_Pt-120to170_MuEnrichedPt5", 5], ["QCD_Pt-15to20_MuEnrichedPt5", 5], ["QCD_Pt-170to300_MuEnrichedPt5", 5], ["QCD_Pt-20to30_MuEnrichedPt5", 5], ["QCD_Pt-300to470_MuEnrichedPt5", 5], ["QCD_Pt-30to50_MuEnrichedPt5", 5], ["QCD_Pt-470to600_MuEnrichedPt5", 5], ["QCD_Pt-50to80_MuEnrichedPt5", 5], ["QCD_Pt-600to800_MuEnrichedPt5", 5], ["QCD_Pt-800to1000_MuEnrichedPt5", 5], ["QCD_Pt-80to120_MuEnrichedPt5", 5],
  #["QCD_Pt_170to250_bcToE", 5], ["QCD_Pt_20to30_bcToE", 5], ["QCD_Pt_250toInf_bcToE", 5], ["QCD_Pt_30to80_bcToE", 5], ["QCD_Pt_80to170_bcToE", 5],
  #["QCD_Pt-120to170_EMEnriched", 5], ["QCD_Pt-15to20_EMEnriched", 5], ["QCD_Pt-20to30_EMEnriched", 5], ["QCD_Pt-300toInf_EMEnriched", 5], ["QCD_Pt-30to50_EMEnriched", 5], ["QCD_Pt-50to80_EMEnriched", 5], ["QCD_Pt-80to120_EMEnriched", 5],
]


List_DiLepValid_Data_2017_DiMuORDiEl = [
  ["DoubleEG", 9, "B"], ["DoubleEG", 21, "C"], ["DoubleEG", 8, "D"], ["DoubleEG", 16, "E"], ["DoubleEG", 21, "F"],
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
]

List_DiLepValid_Data_2017_DiLep = [
  ["DoubleEG", 9, "B"], ["DoubleEG", 21, "C"], ["DoubleEG", 8, "D"], ["DoubleEG", 16, "E"], ["DoubleEG", 21, "F"],
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
  ["MuonEG", 13, "B"], ["MuonEG", 15, "C"], ["MuonEG", 7, "D"], ["MuonEG", 15, "E"], ["MuonEG", 24, "F"],
]

List_DiLepValid_Data_2017_DiElectron = [
#  ["SingleElectron", 9, "B"], ["SingleElectron", 21, "C"], ["SingleElectron", 8, "D"], ["SingleElectron", 16, "E"], ["SingleElectron", 21, "F"],
  ["DoubleEG", 9, "B"], ["DoubleEG", 21, "C"], ["DoubleEG", 8, "D"], ["DoubleEG", 16, "E"], ["DoubleEG", 21, "F"],
]
List_DiLepValid_Data_2017_DiMuon = [
#  ["SingleMuon", 13, "B"], ["SingleMuon", 15, "C"], ["SingleMuon", 7, "D"], ["SingleMuon", 15, "E"], ["SingleMuon", 24, "F"],
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
]
List_DiLepValid_Data_2017_ElectronMuon = [
  ["MuonEG", 13, "B"], ["MuonEG", 15, "C"], ["MuonEG", 7, "D"], ["MuonEG", 15, "E"], ["MuonEG", 24, "F"],
]
List_DiLepValid_Data_2018_DiElectron = [
  ["EGamma", 18, "A"], ["EGamma", 8, "B"], ["EGamma", 8, "C"], ["EGamma", 41, "D"],
]
List_DiLepValid_Data_2018_DiMuon = [
  ["SingleMuon", 18, "A"], ["SingleMuon", 9, "B"], ["SingleMuon", 9, "C"], ["SingleMuon", 39, "D"],
]
List_DiLepValid_BkdMC_2017 = [
  ["DYJets", 21],
  ["TTLL_powheg", 13], ["WZ_pythia", 1], ["WW_pythia", 2], ["ZZ_pythia", 1], ["SingleTop_tW_antitop_NoFullyHad", 2], ["SingleTop_tW_top_NoFullyHad", 2],
]

List_ControlPlots_BkdMC_2017 = [
#  ["DYJets", 21], ["DYJets10to50_MG", 13],
#  ["WZTo2L2Q_amcatnlo", 7], ["WZTo3LNu_powheg", 3], ["ZZTo2L2Q_amcatnlo", 7], ["ZZTo2L2Nu_powheg", 3], ["ZZTo4L_powheg", 3], ["WWTo2L2Nu_powheg", 3], 
#  ["TTLL_powheg", 13], ["SingleTop_tW_antitop_NoFullyHad", 2], ["SingleTop_tW_top_NoFullyHad", 2],
#  ["WZ_pythia", 1], ["WW_pythia", 2], ["ZZ_pythia", 1], 
  ["TTLL_powheg", 13], ["TTLJ_powheg", 7],
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
  #["DYJets", 25], ["TTLL_powheg", 13],
  ["ZZTo4L_powheg", 3], ["ggHToZZTo4L", 2], ["VBF_HToZZTo4L", 1], ["ttHToNonbb", 2], ["VHToNonbb", 1], ["ttZToLLNuNu", 4], ["WWZ", 1], ["WZZ", 1], ["ZZZ", 1],
  ["WZTo3LNu_powheg", 1], ["ttWToLNu", 2], ["WWW", 1], #Only Trilep  
  ["ZGToLLG_01J", 5], ["TTG", 2],
  #["WJets_MG", 7], ["TTLJ_powheg", 7], #Only SS2l
  #["WGToLNuG_01J", 5], ["TG", 1], ["WpWp_EWK", 1], ["WpWp_QCD", 1],#SS2l
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

List_SkimTree_TrigInfo_Data_2017 = [
  ["SingleMuonTrigInfo", 13, "B"], ["SingleMuonTrigInfo", 15, "C"], ["SingleMuonTrigInfo", 7, "D"], ["SingleMuonTrigInfo", 15, "E"], ["SingleMuonTrigInfo", 24, "F"],
  ["SingleElectronTrigInfo", 13, "B"], ["SingleElectronTrigInfo", 15, "C"], ["SingleElectronTrigInfo", 7, "D"], ["SingleElectronTrigInfo", 15, "E"], ["SingleElectronTrigInfo", 24, "F"],
]

List_SkimTree_SS2lOR3l_BkdMC_2017 = [
  ["SingleTop_tW_antitop_NoFullyHad", 5], ["SingleTop_tW_top_NoFullyHad", 5], ["SingleTop_sch_Lep", 5], ["SingleTop_tch_antitop_Incl", 5], ["SingleTop_tch_top_Incl", 5],
]



List_SkimRateCheck_BkdMC_2016 = [
  ["DYJets", 25], ["TTLL_powheg", 13], ["WJets_MG", 13], ["TTLJ_powheg", 25], #Only SS2l
]

List_TrigCheck_BkdMC_2016 = [
  ["WZTo3LNu_powheg", 5], ["ZZTo4L_powheg", 20],
]

List_TrigCheck_BkdMC_2017 = [
  ["TTLL_powheg", 13], #Only SS2l
  ["DYJets", 25],
]

List_MCPUDist17_BkdMC_2017 = [
  ["DYJets_MG_HT-40to70", 10], ["DYJets_MG_HT-70to100", 10], ["DYJets_MG_HT-100to200", 10], ["DYJets_MG_HT-200to400", 10], ["DYJets_MG_HT-400to600", 10], ["DYJets_MG_HT-600to800", 10], ["DYJets_MG_HT-800to1200", 5], ["DYJets_MG_HT-1200to2500", 5], ["DYJets_MG_HT-2500toInf", 5], 
]

List_HNTopV1_SigMC_2017 = [
  ["TT_TTobNMu_SS2L_LO_MN20", 1], ["TT_TTobNMu_SS2L_LO_MN50", 1], ["TT_TTobNMu_SS2L_LO_MN100", 1],
  ["TT_TTobNMu_LepTop3L_LO_MN20", 1], ["TT_TTobNMu_LepTop3L_LO_MN50", 1], ["TT_TTobNMu_LepTop3L_LO_MN100", 1],
  ["TT_TTobNMu_HadTop3L_LO_MN20", 1], ["TT_TTobNMu_HadTop3L_LO_MN50", 1], ["TT_TTobNMu_HadTop3L_LO_MN100", 1],
  ["TT_TTobNMu_4L_LO_MN20", 1], ["TT_TTobNMu_4L_LO_MN50", 1], ["TT_TTobNMu_4L_LO_MN100", 1],
]


List_HNTopV2_SigMC_2017 = [
["TTbarTypeIHeavyN-El_2L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_2L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN75" ,1], 
["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN75" ,1], 
["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN75" ,1], 
["TTbarTypeIHeavyN-El_4L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_4L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN75" ,1],
#["TTbarTypeIHeavyN-El_2L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_2L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN90" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN90" ,1],
#["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-El_LepTop3L_MLM_MN90" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN90" ,1],
#["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-El_HadTop3L_MLM_MN90" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN90" ,1],
#["TTbarTypeIHeavyN-El_4L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_4L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-El_4L_MLM_MN90" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-Mu_4L_MLM_MN90" ,1],
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
#  ["TT_TTobNMu_SS2L_LO_MN20", 1], ["TT_TTobNMu_SS2L_LO_MN50", 1], ["TT_TTobNMu_SS2L_LO_MN100", 1],
  ["TTbarTypeIHeavyN-El_2L_MLM_MN100",1], ["TTbarTypeIHeavyN-El_2L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-El_2L_MLM_MN90" ,1],
#  ["TTbarTypeIHeavyN-Mu_2L_MLM_MN100",1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN20" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN30" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN40" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN50" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN60" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN70" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN75" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN85" ,1], ["TTbarTypeIHeavyN-Mu_2L_MLM_MN90" ,1]
]

List_IDOptimization_BkdMC_2016 = [
  ["TTLJ_powheg", 25],
]

List_IDOptimization_BkdMC_2017 = [
  ["ZGToLLG_01J", 5],
  ["TTLL_powheg", 13],
  ["TTLJ_powheg", 15],
#  ["TTLL_powheg_TrigInfo", 13],
#  ["TTLJ_powheg_TrigInfo", 15],
#  ["DYJets", 25],
#  ["WJets_MG", 7],
#  ["DYJets_TrigInfo", 10],
]

List_IDOptimization_BkdMC_2018 = [
  ["TTLJ_powheg", 25],
]


List_MeasTrigEff_BkdMC_2017 = [["TTLL_powheg_TrigInfo", 13], ["DYJets_TrigInfo", 10]]
#List_MeasTrigEff_BkdMC_2017 = [["TTLL_powheg", 13], ["DYJets",25]]
List_MeasTrigEff_Data_2017_DiMu = [
  ["DoubleMuon", 13, "B"], ["DoubleMuon", 15, "C"], ["DoubleMuon", 7, "D"], ["DoubleMuon", 15, "E"], ["DoubleMuon", 24, "F"],
]
List_MeasTrigEff_Data_2017_ElMu = [
  ["MuonEG", 13, "B"], ["MuonEG", 15, "C"], ["MuonEG", 7, "D"], ["MuonEG", 15, "E"], ["MuonEG", 24, "F"],
#  ["SingleMuonTrigInfo", 13, "B"], ["SingleMuonTrigInfo", 15, "C"], ["SingleMuonTrigInfo", 7, "D"], ["SingleMuonTrigInfo", 15, "E"], ["SingleMuonTrigInfo", 24, "F"],
#  ["SingleElectronTrigInfo", 13, "B"], ["SingleElectronTrigInfo", 15, "C"], ["SingleElectronTrigInfo", 7, "D"], ["SingleElectronTrigInfo", 15, "E"], ["SingleElectronTrigInfo", 24, "F"],
]






List_Debug_SigMC_2016 = [["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1]]
List_Debug_BkdMC_2016 = [
  ["TTLJ_powheg", 1]
  #["ZGTo2LG", 1]
]
List_Debug_Data_2016 = [["MuonEG", 1, "B_ver2"]]
List_Debug_SigMC_2017 = [
  #["TT_TTobNMu_SS2L_LO_MN50", 1]
  #["TT_TTobNMu_LepTop3L_LO_MN20", 1]
  #["TT_TTobNMu_4L_LO_MN20", 1]
#  ["TTbarTypeIHeavyN-El_2L_MLM_MN20", 1],
#  ["TTbarTypeIHeavyN-Mu_2L_MLM_MN100", 1],
#  ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN20", 1],
#  ["TTbarTypeIHeavyN-Mu_HadTop3L_MLM_MN100", 1],
#  ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN20", 1],
#  ["TTbarTypeIHeavyN-Mu_LepTop3L_MLM_MN100", 1],
#  ["TTbarTypeIHeavyN-Mu_4L_MLM_MN20", 1],
#  ["TTbarTypeIHeavyN-Mu_4L_MLM_MN100", 1],
]
List_Debug_BkdMC_2017 = [["DYJets", 1]]
#List_Debug_BkdMC_2017 = [["ZGToLLG_01J", 1]]
#List_Debug_BkdMC_2017 = [["TTLL_powheg", 1]]
List_Debug_Data_2017 = [["DoubleEG", 1, "B"]]
#List_Debug_Data_2017 = [["DoubleMuon", 1, "B"]]
List_Debug_SigMC_2018 = [["TestFiltered_TTToHcToWA_AToMuMu_MHc160_MA75", 1]]
List_Debug_BkdMC_2018 = [["ZZTo4L_powheg", 1]]
List_Debug_Data_2018 = [["SingleMuon", 1, "B"]]



Dict_JobConfig = {
 "MC_NonFk3l_2017": List_MC_NonFk3l_2017,
 "FakeRateMeas_Data_2017_El": List_FakeRateMeas_Data_2017_El,
 "FakeRateMeas_Data_2017_Mu": List_FakeRateMeas_Data_2017_Mu,
 "FakeRateMeas_BkdMC_2017": List_FakeRateMeas_BkdMC_2017,
 "MeasTrigEff_Data_2017_ElectronMuon": List_MeasTrigEff_Data_2017_ElMu,
 "MeasTrigEff_Data_2017_DiMuon": List_MeasTrigEff_Data_2017_DiMu,
 "MeasTrigEff_BkdMC_2017": List_MeasTrigEff_BkdMC_2017,
 "IDOptimization_BkdMC_2016": List_IDOptimization_BkdMC_2016,
 "IDOptimization_BkdMC_2017": List_IDOptimization_BkdMC_2017,
 "IDOptimization_BkdMC_2018": List_IDOptimization_BkdMC_2018,
 "IDOptimization_SigMC_2017": List_GenMatchingValid_SigMC_2017,
 "GenMatchingValid_BkdMC_2016": List_GenMatchingValid_BkdMC_2016,
 "GenMatchingValid_BkdMC_2017": List_GenMatchingValid_BkdMC_2017,
 "GenMatchingValid_BkdMC_2018": List_GenMatchingValid_BkdMC_2018,
 "GenMatchingValid_SigMC_2017": List_GenMatchingValid_SigMC_2017,
 "TrigCheck_Data_2017_DiMuon": List_TestRun_Data_2017_TriMu,
 "TrigCheck_BkdMC_2016": List_TrigCheck_BkdMC_2016,
 "TrigCheck_BkdMC_2017": List_TrigCheck_BkdMC_2017,
 "SkimRateCheck_BkdMC_2016": List_SkimRateCheck_BkdMC_2016,
 "SkimRateCheck_Data_2016_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2016,
 "SkimRateCheck_Data_2017_TrigInfo": List_SkimTree_TrigInfo_Data_2017,
 "SkimTree_SS2lOR3l_BkdMC_2016": List_HNTopFeas_BkdMC_2016,
 "SkimTree_SS2lOR3l_BkdMC_2017": List_SkimTree_SS2lOR3l_BkdMC_2017,
 "SkimTree_SS2lOR3l_BkdMC_2018": List_HNTopFeas_BkdMC_2018,
 "SkimTree_SS2lOR3l_SigMC_2017": List_HNTopV2_SigMC_2017,
 "SkimTree_SS2lOR3l_Data_2016_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2016,
 "SkimTree_SS2lOR3l_Data_2017_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2017,
 "SkimTree_SS2lOR3l_Data_2018_SS2lOR3l": List_SkimTree_SS2lOR3l_Data_2018,
 "SyncYield_Data_2016_ElDiMu": List_TestRun_Data_2016_ElDiMu,
 "SyncYield_Data_2016_TriMu": List_TestRun_Data_2016_TriMu,
 "MCPUDist17_BkdMC_2017" : List_MCPUDist17_BkdMC_2017,
 #"MeasCFlipRate_BkdMC_2017": List_MC_Fk3l_2017,
 "MeasCFlipRate_BkdMC_2017": List_MC_CF2l_2017,
 #"MeasCFlipRate_BkdMC_2017": List_MC_CFMeas_2017,
 "MeasCFlipRate_Data_2017_All2l": List_Data_DiEl_2017,
 "ControlPlots_BkdMC_2017": List_MC_Cv2l_2017,
 #"ControlPlots_BkdMC_2017": List_MC_CF2l_2017,
 #"ControlPlots_BkdMC_2017": List_MC_CF2l_2017,
 #"ControlPlots_BkdMC_2017": List_MC_Pr2l_2017,
 #"ControlPlots_BkdMC_2017": List_MC_Cv3l_2017,
 #"ControlPlots_BkdMC_2017": List_MC_NonFk3l_2017,
 "ControlPlots_SigMC_2017": List_HNTopV2_SigMC_2017,
 "ControlPlots_Data_2017_OS2l": List_DiLepValid_Data_2017_DiMuORDiEl,
 "ControlPlots_Data_2017_SS2l": List_Data_SF2L_2017,
 "ControlPlots_Data_2017_TriLep": List_Data_2LPD_2017,
 "KinVarSearch_BkdMC_2017": List_HNTopFeas_BkdMC_2017,
 "KinVarSearch_SigMC_2017": List_HNTopV2_SigMC_2017,
 "KinVarPlot_Conv_2017": List_MC_Cv2l_2017,
 "KinVarPlot_BkdMC_2017": List_HNTopFeas_BkdMC_2017,
 "KinVarPlot_SigMC_2017": List_HNTopV2_SigMC_2017,
 "KinVarPlot_Data_2017_SS2l": List_Data_SF2L_2017,
 "KinVarPlot_Flip_2017_SS2l": List_Data_DiEl_2017,
 "HNTopFeas_BkdMC_2016": List_HNTopFeas_BkdMC_2016,
 "HNTopFeas_BkdMC_2017": List_HNTopFeas_BkdMC_2017,
 "HNTopFeas_BkdMC_2018": List_HNTopFeas_BkdMC_2018,
 "HNTopFeas_SigMC_2017": List_HNTopV2_SigMC_2017,
 "TestRun_Data_2016_TriMu": List_TestRun_Data_2016_TriMu,
 "TestRun_Data_2017_TriMu": List_TestRun_Data_2017_TriMu,
 "TestRun_Data_2018_TriMu": List_TestRun_Data_2018_TriMu,
 "TestRun_Data_2016_ElDiMu": List_TestRun_Data_2016_ElDiMu,
 "TestRun_Data_2017_ElDiMu": List_TestRun_Data_2017_ElDiMu,
 "TestRun_Data_2018_ElDiMu": List_TestRun_Data_2018_ElDiMu,
 "DiLepValid_Data_2017_DiElectron": List_DiLepValid_Data_2017_DiElectron,
 "DiLepValid_Data_2017_DiMuon": List_DiLepValid_Data_2017_DiMuon,
 "DiLepValid_Data_2017_ElectronMuon": List_DiLepValid_Data_2017_ElectronMuon,
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
