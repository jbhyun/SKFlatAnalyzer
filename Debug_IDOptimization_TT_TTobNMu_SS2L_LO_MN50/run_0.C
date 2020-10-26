R__LOAD_LIBRARY(/data8/Users/jbhyun/cmssw/SKFlatAnalyzer/SKFlatAnalyzer/lib/libAnalyzers.so)
R__LOAD_LIBRARY(/data8/Users/jbhyun/cmssw/SKFlatAnalyzer/SKFlatAnalyzer/lib/libDataFormats.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/lib/libLHAPDF.so)


void run_0(){

  IDOptimization m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 10000;
  m.MaxEvent = 100000;
  m.MCSample = "TT_TTobNMu_SS2L_LO_MN50";
  //m.MCSample = "Trig_TT_TTobNMu_SS2L_LO_MN50";
  m.IsDATA = false;
  m.xsec = 1;
  m.sumW = 192000;
  m.IsFastSim = false;
  m.DataYear = 2017;
    //"MuID",
    //"MuTrigEffect", 
  m.Userflags = {
    //"MuMu", "TrigEffCheck", 
    "MuMu", "SelEffCheck", 
  };
  //m.AddFile("/data8/Users/jbhyun/cmssw/SKFlatMaker/CMSSW_10_2_18/src/SKFlatMaker/SKFlatMaker/script/TAMSA/SKFlatNtuple_2017_MC_0.root");
  //m.AddFile("/data8/Users/jbhyun/cmssw/SKFlatMaker/CMSSW_10_2_18/src/SKFlatMaker/SKFlatMaker/script/TAMSA/SKFlatNtuple_2017_MC_1.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_24.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_55.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_75.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_78.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_79.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_36.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_37.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_44.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_47.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_56.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_91.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_45.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_52.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_54.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_6.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_69.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_76.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_87.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_11.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_13.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_16.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_3.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_35.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_39.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_40.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_41.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_53.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_68.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_80.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_82.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_92.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_93.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_94.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_96.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_1.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_14.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_25.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_28.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_33.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_4.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_42.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_58.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_65.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_66.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_8.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_85.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_88.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_9.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_15.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_17.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_18.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_21.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_22.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_49.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_51.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_60.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_62.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_71.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_72.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_81.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_95.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_10.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_12.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_2.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_20.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_23.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_26.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_27.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_29.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_30.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_31.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_32.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_34.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_38.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_43.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_46.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_48.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_5.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_50.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_57.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_59.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_61.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_63.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_64.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_67.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_7.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_70.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_73.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_74.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_77.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_83.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_84.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_86.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_89.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_90.root");
//  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/preliminary_TTbarTypeIHeavuN-Mu/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/SKFlat_Run2Legacy_v4/200718_091024/0000/SKFlatNtuple_2017_MC_19.root");
  m.SetOutfilePath("hists.root");
  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}
