#include "BTagEff_MC.h"

void BTagEff_MC::initializeAnalyzer(){

  MCEffMeas=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("MCEffMeas")) MCEffMeas=true;
    if(Userflags.at(i).Contains("SystRun")) SystRun=true; 
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void BTagEff_MC::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  if(!PassMETFilter()) return;
  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);

  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "POGIDTIsoVVL", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "passMVAID_iso_WP90", 10., 2.5);


  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  std::vector<Jet> jetNoVetoColl  = GetJets("tight", 20., 2.4);
  //std::vector<Jet> bjetNoVetoColl = SelBJets(jetNoVetoColl, param_jets);
  std::vector<Jet> jetColl  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl, muonLooseColl, 0.4);
  std::vector<Jet> jetNoPrColl;
  //std::vector<Jet> bjetColl = SelBJets(jetColl, param_jets);


  std::vector<Gen> truthColl=GetGens();
  std::vector<int> IdxColl_PrLep;
  for(unsigned int it=2; it<truthColl.size(); it++){
    int LepType = GetLeptonType_JH(it, truthColl);
    if(LepType>0 && LepType<4) IdxColl_PrLep.push_back(it);
  }
  for(unsigned int it_j=0; it_j<jetNoVetoColl.size(); it_j++){
    bool ContainsPrLep=false;
    for(unsigned int it_l=0; it_l<IdxColl_PrLep.size(); it_l++){
      FillHist("dRlj", jetNoVetoColl.at(it_j).DeltaR(truthColl.at(IdxColl_PrLep.at(it_l))), weight, 50, 0., 5.);
      FillHist("Etal", truthColl.at(IdxColl_PrLep.at(it_l)).Eta(), weight, 20, -5., 5.);
      FillHist("Ptl", truthColl.at(IdxColl_PrLep.at(it_l)).Pt(), weight, 100, 0., 500.);
      if(jetNoVetoColl.at(it_j).DeltaR(truthColl.at(IdxColl_PrLep.at(it_l)))<0.4){ ContainsPrLep=true; break; }
    }
    if(!ContainsPrLep) jetNoPrColl.push_back(jetNoVetoColl.at(it_j));
  }
  if(jetNoPrColl.size()==jetNoVetoColl.size()) FillHist("Check", 0., weight, 5, 0., 5.);
  else FillHist("Check", 1., weight, 5, 0., 5.);
  if(IdxColl_PrLep.size()!=0) FillHist("Check", 2., weight, 5, 0., 5);


  bool EventCand = true;
  float w_gen = 1., w_filter = 1., w_topptrw = 1., w_lumi = 1., w_PU = 1., w_prefire = 1., sf_trig = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
  if((!IsDATA) and EventCand){
    w_gen     = ev.MCweight();
    //w_filter  = GetGenFilterEffCorr();
    //w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    w_lumi    = weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    w_PU      = GetPileUpWeight(nPileUp, 0);
    w_prefire = GetPrefireWeight(0);
    //sf_muid   = GetMuonSF(muonTightColl, "POGTID_genTrk", "ID");
    //sf_muiso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    //sf_elreco = GetElectronSF(electronTightColl, "", "Reco");
    //sf_elid   = GetElectronSF(electronTightColl, "POGMVAIsoWP90", "ID");
    //sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    //sf_trig   = mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig, "");
    //cout<<"w_gen:"<<w_gen<<" w_lumi:"<<w_lumi<<" w_PU:"<<w_PU<<" w_prefire:"<<w_prefire<<" sf_trig:"<<sf_trig<<endl;
    //cout<<"sf_mutk"<<sf_mutk<<" sf_muid:"<<sf_muid<<" sf_muiso:"<<sf_muiso<<" sf_elreco:"<<sf_elreco<<" sf_elid:"<<sf_elid<<" sf_btag:"<<sf_btag<<endl;
  }
  weight *= w_gen * w_filter * w_topptrw * w_lumi * w_PU * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

 
  if(MCEffMeas){
    MeasureMCBTagEfficiency(jetNoVetoColl, param_jets, weight, "_Incl");
    MeasureMCBTagEfficiency(jetColl, param_jets, weight, "_LepVeto");
    MeasureMCBTagEfficiency(jetNoPrColl, param_jets, weight, "_PrLVeto");
  }

}



void BTagEff_MC::MeasureMCBTagEfficiency(std::vector<Jet>& JetColl, JetTagging::Parameters jtp, float weight, TString Label)

{

  const int NPtBinEdges=10, NfEtaBinEdges=5;
  double fEtaBinEdges[NfEtaBinEdges]={0., 0.6, 1.2, 1.8, 2.4};
  double PtBinEdges[NPtBinEdges]={20., 30., 50., 70., 100., 140., 200., 300., 600., 1000.};

  for(unsigned int it_j=0; it_j<JetColl.size(); it_j++){
    double PTj   = JetColl.at(it_j).Pt();
    double fEtaj = fabs(JetColl.at(it_j).Eta());
    bool   IsTagged = JetColl.at(it_j).GetTaggerResult(jtp.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtp.j_Tagger, jtp.j_WP);
    int    JetFlav  = abs(JetColl.at(it_j).hadronFlavour());
    if(JetFlav==5){
      FillHist("Nb_PtEta"+Label, PTj, fEtaj, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
      if(IsTagged) FillHist("NbTag_PtEta"+Label, PTj, fEtaj, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
    }
    else if(JetFlav==4){
      FillHist("Nc_PtEta"+Label, PTj, fEtaj, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
      if(IsTagged) FillHist("NcTag_PtEta"+Label, PTj, fEtaj, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
    }
    else if(JetFlav==0){
      FillHist("Nl_PtEta"+Label, PTj, fEtaj, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
      if(IsTagged) FillHist("NlTag_PtEta"+Label, PTj, fEtaj, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
    }
    else FillHist("NError"+Label, 0., weight, 1, 0., 1.);
  }

}





void BTagEff_MC::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

BTagEff_MC::BTagEff_MC(){

}

BTagEff_MC::~BTagEff_MC(){

}


