import os, subprocess, re 
from os import path


if not "SKFlat_WD" in os.environ:
  print("Set up SKFlat environment."); exit(); 

#[PD(Or PrivateMC dir name), alias, xsec(fb)]
Arr_Sample = [["preliminary_TTbarTypeIHeavyN-Mu_LepTop3L_LO_MN20" , "TT_TTobNMu_LepTop3L_LO_MN20" , 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_LepTop3L_LO_MN50" , "TT_TTobNMu_LepTop3L_LO_MN50" , 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_LepTop3L_LO_MN100", "TT_TTobNMu_LepTop3L_LO_MN100", 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_HadTop3L_LO_MN20" , "TT_TTobNMu_HadTop3L_LO_MN20" , 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_HadTop3L_LO_MN50" , "TT_TTobNMu_HadTop3L_LO_MN50" , 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_HadTop3L_LO_MN100", "TT_TTobNMu_HadTop3L_LO_MN100", 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_4L_LO_MN20" , "TT_TTobNMu_4L_LO_MN20" , 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_4L_LO_MN50" , "TT_TTobNMu_4L_LO_MN50" , 1.],
              ["preliminary_TTbarTypeIHeavyN-Mu_4L_LO_MN100", "TT_TTobNMu_4L_LO_MN100", 1.]]


VerProc="Run2Legacy_v4"
Year="2017"
DataType="PrivateMC"
Category="preliminary_TTbarTypeIHeavuN-Mu" #"HctoWA"
#SamplePath="/data9/Users/jbhyun/SKFlat/"+VerProc+"/"+Year+"/"+DataType+"/"+Category+"/" #HctoWA Test
#SamplePath="/gv0/Users/jbhyun/SKFlat/TrigStudy/"+Year+"/" #TrigStudy
SamplePath="/gv0/DATA/SKFlat/"+VerProc+"/"+Year+"/"+DataType+"/"+Category+"/" #Normal samples
OverWrite = True
Verbose = True

WorkingPath=os.getenv('SKFlat_WD')
ListPath=WorkingPath+"/data/"+VerProc+"/"+Year+"/Sample/"
CalcCodePath= WorkingPath+"/script/MakeSampleInfo/"

for it_proc in Arr_Sample:
  InputDirName = it_proc[0]
  Alias        = it_proc[1]
  Xsec         = it_proc[2]
  
  Path_Smplinfo = ListPath+"CommonSampleInfo/"+Alias+".txt"
  Path_SmplPath = ListPath+"ForSNU/"+Alias+".txt"

  if path.exists(Path_SmplPath):
    if not OverWrite:
      print("[Error] Skip due to existing path: "+Path_SmplPath); continue;
    else:
      os.system("rm "+Path_SmplPath)
  if path.exists(Path_Smplinfo):
    if not OverWrite:
      print("[Error] Skip due to existing path: "+Path_Smplinfo); continue;
    else:
      os.system("rm "+Path_Smplinfo)


  os.system("touch "+Path_SmplPath)
  os.system("find "+SamplePath+InputDirName+" -name *root >> "+Path_SmplPath)
  os.system("touch "+Path_Smplinfo)
  os.system("echo \"# alias PD xsec nmc sumw\" >> "+Path_Smplinfo)
  CountInfo = subprocess.check_output("root -l -b -q \""+CalcCodePath+"LumiWCalc.C(\\\""+Path_SmplPath+"\\\")\"", shell=True);
  Nevent = CountInfo.split()[2]
  SumW   = CountInfo.split()[3]
  os.system("echo \""+Alias+"\t"+InputDirName+"\t"+str(Xsec)+"\t"+Nevent+"\t"+SumW+"\" >> "+Path_Smplinfo)

  if Verbose:
    print("Processed "+Alias+" ("+VerProc+", "+Year+", "+DataType+").")

