import os, subprocess, re 
from os import path


if not "SKFlat_WD" in os.environ:
  print("Set up SKFlat environment."); exit(); 

#[PD(Or PrivateMC dir name), alias, xsec(fb)]
#Arr_MCSample = [#16
#  "DYJets", "TTLL_powheg", "ZZTo4L_powheg", "ggHToZZTo4L", "VBF_HToZZTo4L", "ttHToNonbb", "ttZToLLNuNu", "WWZ", "WZZ", "ZZZ",
#  "WZTo3LNu_powheg", "ttWToLNu", "WWW",  #Only Trilep  
#  "WJets_MG", "TTLJ_powheg", #Only SS2l
#  "WGToLNuG", "ZGTo2LG", "TTG", "TG", "VHToNonbb", "WpWp_EWK", "WpWp_QCD", 
#]
#Arr_MCSample = [#17
#  "DYJets", "TTLL_powheg", 
#  "ZZTo4L_powheg", "ggHToZZTo4L", "VBF_HToZZTo4L", "ttHToNonbb", "ttZToLLNuNu", "WWZ", "WZZ", "ZZZ", 
#  "WZTo3LNu_powheg", "ttWToLNu", "WWW", #Only Trilep  
#  "WJets_MG", "TTLJ_powheg", #Only SS2l
#  "WGToLNuG_01J", "ZGToLLG_01J", "TTG", "TG", "VHToNonbb", "WpWp_EWK", "WpWp_QCD",
#]
#Arr_MCSample =["TT_TTobNMu_SS2L_LO_MN20", "TT_TTobNMu_SS2L_LO_MN50", "TT_TTobNMu_SS2L_LO_MN100",
#               "TT_TTobNMu_LepTop3L_LO_MN20", "TT_TTobNMu_LepTop3L_LO_MN50", "TT_TTobNMu_LepTop3L_LO_MN100",
#               "TT_TTobNMu_HadTop3L_LO_MN20", "TT_TTobNMu_HadTop3L_LO_MN50", "TT_TTobNMu_HadTop3L_LO_MN100",
#               "TT_TTobNMu_4L_LO_MN20", "TT_TTobNMu_4L_LO_MN50", "TT_TTobNMu_4L_LO_MN100"]
Arr_MCSample =["QCD_Pt-1000toInf_MuEnrichedPt5", "QCD_Pt-1000toInf_MuEnrichedPt5", "QCD_Pt-1000toInf_MuEnrichedPt5", "QCD_Pt-120to170_EMEnriched", "QCD_Pt-120to170_EMEnriched", "QCD_Pt-120to170_EMEnriched", "QCD_Pt-120to170_MuEnrichedPt5", "QCD_Pt-120to170_MuEnrichedPt5", "QCD_Pt-120to170_MuEnrichedPt5", "QCD_Pt-15to20_EMEnriched", "QCD_Pt-15to20_EMEnriched", "QCD_Pt-15to20_EMEnriched", "QCD_Pt-15to20_MuEnrichedPt5", "QCD_Pt-15to20_MuEnrichedPt5", "QCD_Pt-15to20_MuEnrichedPt5", "QCD_Pt-170to300_MuEnrichedPt5", "QCD_Pt-170to300_MuEnrichedPt5", "QCD_Pt-170to300_MuEnrichedPt5", "QCD_Pt-20to30_EMEnriched", "QCD_Pt-20to30_EMEnriched", "QCD_Pt-20to30_EMEnriched", "QCD_Pt-20to30_MuEnrichedPt5", "QCD_Pt-20to30_MuEnrichedPt5", "QCD_Pt-20to30_MuEnrichedPt5", "QCD_Pt-300to470_MuEnrichedPt5", "QCD_Pt-300to470_MuEnrichedPt5", "QCD_Pt-300to470_MuEnrichedPt5", "QCD_Pt-300toInf_EMEnriched", "QCD_Pt-300toInf_EMEnriched", "QCD_Pt-300toInf_EMEnriched", "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf", "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf", "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf", "QCD_Pt-30to50_EMEnriched", "QCD_Pt-30to50_EMEnriched", "QCD_Pt-30to50_EMEnriched", "QCD_Pt-30to50_MuEnrichedPt5", "QCD_Pt-30to50_MuEnrichedPt5", "QCD_Pt-30to50_MuEnrichedPt5", "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf", "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf", "QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf", "QCD_Pt-470to600_MuEnrichedPt5", "QCD_Pt-470to600_MuEnrichedPt5", "QCD_Pt-470to600_MuEnrichedPt5", "QCD_Pt-50to80_EMEnriched", "QCD_Pt-50to80_EMEnriched", "QCD_Pt-50to80_EMEnriched", "QCD_Pt-50to80_MuEnrichedPt5", "QCD_Pt-50to80_MuEnrichedPt5", "QCD_Pt-50to80_MuEnrichedPt5", "QCD_Pt-600to800_MuEnrichedPt5", "QCD_Pt-600to800_MuEnrichedPt5", "QCD_Pt-600to800_MuEnrichedPt5", "QCD_Pt-800to1000_MuEnrichedPt5", "QCD_Pt-800to1000_MuEnrichedPt5", "QCD_Pt-800to1000_MuEnrichedPt5", "QCD_Pt-80to120_EMEnriched", "QCD_Pt-80to120_EMEnriched", "QCD_Pt-80to120_EMEnriched", "QCD_Pt-80to120_MuEnrichedPt5", "QCD_Pt-80to120_MuEnrichedPt5", "QCD_Pt-80to120_MuEnrichedPt5", "QCD_Pt_170to250_bcToE", "QCD_Pt_170to250_bcToE", "QCD_Pt_170to250_bcToE", "QCD_Pt_20to30_bcToE", "QCD_Pt_20to30_bcToE", "QCD_Pt_20to30_bcToE", "QCD_Pt_250toInf_bcToE", "QCD_Pt_250toInf_bcToE", "QCD_Pt_250toInf_bcToE", "QCD_Pt_30to80_bcToE", "QCD_Pt_30to80_bcToE", "QCD_Pt_30to80_bcToE", "QCD_Pt_80to170_bcToE", "QCD_Pt_80to170_bcToE", "QCD_Pt_80to170_bcToE"]


#Arr_MCSample = [#18
#  "DYJets", "TTLL_powheg", 
#  "ZZTo4L_powheg", "ggHToZZTo4L", "VBF_HToZZTo4L", "ttHToNonbb", "ttZToLLNuNu", "WWZ", "WZZ", "ZZZ", 
#  "WZTo3LNu_powheg", "ttWToLNu", "WWW",  #Only Trilep  
#  "WJets_MG", "TTLJ_powheg",  #Only SS2l
#  "WGToLNuG_01J", "ZGToLLG_01J", "TTG", "TG", "VHToNonbb", "WpWp_EWK", "WpWp_QCD", 
#]
Arr_DataSample = [#16
  ["DoubleMuon", "B_ver2"], ["DoubleMuon", "C"], ["DoubleMuon", "D"], ["DoubleMuon", "E"], ["DoubleMuon", "F"], ["DoubleMuon", "G"], ["DoubleMuon", "H"],
  ["MuonEG", "B_ver2"], ["MuonEG", "C"], ["MuonEG", "D"], ["MuonEG", "E"], ["MuonEG", "F"], ["MuonEG", "G"], ["MuonEG", "H"],
  ["DoubleEG", "B_ver2"], ["DoubleEG", "C"], ["DoubleEG", "D"], ["DoubleEG", "E"], ["DoubleEG", "F"], ["DoubleEG", "G"], ["DoubleEG", "H"],
]
#Arr_DataSample = [#17
#  ["DoubleMuon", "B"], ["DoubleMuon", "C"], ["DoubleMuon", "D"], ["DoubleMuon", "E"], ["DoubleMuon", "F"],
#  ["MuonEG", "B"], ["MuonEG", "C"], ["MuonEG", "D"], ["MuonEG", "E"], ["MuonEG", "F"],
#  ["DoubleEG", "B"], ["DoubleEG", "C"], ["DoubleEG", "D"], ["DoubleEG", "E"], ["DoubleEG", "F"],
#]
#Arr_DataSample = [#18
#  ["DoubleMuon", "A"], ["DoubleMuon", "B"], ["DoubleMuon", "C"], ["DoubleMuon", "D"],
#  ["MuonEG", "A"], ["MuonEG", "B"], ["MuonEG", "C"], ["MuonEG", "D"],
#  ["EGamma", "A"], ["EGamma", "B"], ["EGamma", "C"], ["EGamma", "D"],
#]



VerProc="Run2Legacy_v4"
Year="2017"
DataType="MC" #MC/DATA
#DataType="DATA" #MC/DATA
SkimName="SS2lOR3l"
SamplePath="/gv0/DATA/SKFlat/"+VerProc+"/"+Year+"/"+DataType+"_SkimTree_"+SkimName+"/"
OverWrite = True
Verbose = True

WorkingPath=os.getenv('SKFlat_WD')
ListPath=WorkingPath+"/data/"+VerProc+"/"+Year+"/Sample/"

if DataType == "MC":
  for Alias in Arr_MCSample:
    Path_Smplinfo = ListPath+"CommonSampleInfo/"+Alias+".txt"
    Path_SmplPath = ListPath+"ForSNU/SkimTree_"+SkimName+"_"+Alias+".txt"
  
    if not path.exists(Path_Smplinfo):
      print("[Error] Skip due to existing path: "+Path_Smplinfo); continue;
    if path.exists(Path_SmplPath):
      if not OverWrite:
        print("[Error] Skip due to existing path: "+Path_SmplPath); continue;
      else:
        os.system("rm "+Path_SmplPath)
  
    SmplinfoLine = subprocess.check_output("sed -n -e '2p' "+Path_Smplinfo, shell=True);
    SmplinfoArray = SmplinfoLine.split();
    if SmplinfoArray[0] != Alias:
      print("[Error] Sample info not found in the expected location: "+Alias); continue;
  
    DatasetName = SmplinfoArray[1]
    os.system("touch "+Path_SmplPath)
    os.system("find "+SamplePath+DatasetName+" -name *root >> "+Path_SmplPath)
  
  
    if Verbose:
      print("Processed "+Alias+" ("+VerProc+", "+Year+", MC).")

if DataType == "DATA":
  for it_proc in Arr_DataSample:
    Stream = it_proc[0]
    Period = it_proc[1]
    DatasetName = Stream+"_"+Period
    Path_SmplPath = ListPath+"ForSNU/SkimTree_"+SkimName+"_"+DatasetName+".txt"

    if not path.exists(SamplePath+Stream+"/period"+Period):
      print("[Error] Skip due to non-existing path for "+DatasetName); continue;
    if path.exists(Path_SmplPath):
      if not OverWrite:
        print("[Error] Skip due to existing path: "+Path_SmplPath); continue;
      else:
        os.system("rm "+Path_SmplPath)

    os.system("touch "+Path_SmplPath)
    os.system("find "+SamplePath+Stream+"/period"+Period+" -name *root >> "+Path_SmplPath)
  

    if Verbose:
      print("Processed "+DatasetName+" ("+VerProc+", "+Year+", Data).")
