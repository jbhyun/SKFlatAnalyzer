R__LOAD_LIBRARY(/data8/Users/jbhyun/cmssw/SKFlatAnalyzer/SKFlatAnalyzer/lib/libAnalyzers.so)
R__LOAD_LIBRARY(/data8/Users/jbhyun/cmssw/SKFlatAnalyzer/SKFlatAnalyzer/lib/libDataFormats.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/lib/libLHAPDF.so)


void run_0(){

  HNTopFeas m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 10000;
  m.MaxEvent = 100000;
  m.MCSample = "TTLJ_powheg";
  m.IsDATA = false;
  m.xsec = 365.34;
  m.sumW = 43379133;
  m.IsFastSim = false;
  m.DataYear = 2017;
  m.Userflags = {
    "SS2l",
  };
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_9.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_8.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_1.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_2.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_12.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_10.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_0.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_5.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_7.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_6.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_13.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_11.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_4.root");
  m.AddFile("/gv0/DATA/SKFlat/Run2Legacy_v4/2017/MC_SkimTree_SS2lOR3l/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2020_07_09_130236/SKFlatNtuple_2017_MC_3.root");
  m.SetOutfilePath("hists.root");
  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}
