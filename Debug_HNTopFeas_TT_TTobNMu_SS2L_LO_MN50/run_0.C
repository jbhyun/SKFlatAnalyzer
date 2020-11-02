R__LOAD_LIBRARY(/data8/Users/jbhyun/cmssw/SKFlatAnalyzer/SKFlatAnalyzer/lib/libAnalyzers.so)
R__LOAD_LIBRARY(/data8/Users/jbhyun/cmssw/SKFlatAnalyzer/SKFlatAnalyzer/lib/libDataFormats.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/lib/libLHAPDF.so)


void run_0(){

  HNTopFeas m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 10000;
  m.MaxEvent = 10000;
  m.MCSample = "TT_TTobNMu_SS2L_LO_MN50";
  m.IsDATA = false;
  m.xsec = 1;
  m.sumW = 192000;
  m.IsFastSim = false;
  m.DataYear = 2017;
  m.Userflags = {
    "SS2l",
  };
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/2020_10_21_200506/SKFlatNtuple_2017_MC_1.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/preliminary_TTbarTypeIHeavyN-Mu_SS2L_LO_MN50/2020_10_21_200506/SKFlatNtuple_2017_MC_0.root");
  m.SetOutfilePath("hists.root");
  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}
