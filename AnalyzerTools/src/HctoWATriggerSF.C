#include "MCCorrection.h"

float MCCorrection::GetTriggerSF(vector<Electron>& EleColl, vector<Muon>& MuColl, TString SFKey, TString Option){

  if(IsDATA) return 1.;

  TString NominalOpt=Option; NominalOpt.ReplaceAll("Syst","");
  bool SystRun=Option.Contains("Syst");
  float SystDir=0., RelSystData=0., RelSystMC=0., TotRelSyst=0.;

  float TriggerEff_Data = TriggerEfficiency(EleColl, MuColl, SFKey, true,  NominalOpt);
  float TriggerEff_MC   = TriggerEfficiency(EleColl, MuColl, SFKey, false, NominalOpt);
  if(SystRun){
    float TriggerEff_Data_syst = TriggerEfficiency(EleColl, MuColl, SFKey, true,  Option);
    float TriggerEff_MC_syst   = TriggerEfficiency(EleColl, MuColl, SFKey, false, Option);
    RelSystData = TriggerEff_Data!=0.? (TriggerEff_Data_syst-TriggerEff_Data)/TriggerEff_Data:0.;
    RelSystMC   = TriggerEff_MC  !=0.? (TriggerEff_MC_syst  -TriggerEff_MC  )/TriggerEff_MC  :0.;
    TotRelSyst  = sqrt(pow(RelSystData,2.)+pow(RelSystMC,2.));
    if     (Option.Contains("Up"))   SystDir =  1.;
    else if(Option.Contains("Down")) SystDir = -1.;
  }
  
  float TriggerScaleFactor = TriggerEff_MC!=0.? TriggerEff_Data/TriggerEff_MC:0.;
   if(TriggerScaleFactor<0) TriggerScaleFactor=0.;
   if(SystRun) TriggerScaleFactor *= (1.+SystDir*TotRelSyst);

  return TriggerScaleFactor;

}


float MCCorrection::TriggerEfficiency(vector<Electron>& EleColl, vector<Muon>& MuColl, TString SFKey, bool ReturnDataEff, TString Option){
  //DataorMC : T: Return DataEff, F: Return MCEff

  if(IsDATA) return 1.;

  TString StrMCorData = ReturnDataEff? "DATA":"MC";
  int SystDir=0;
  if(Option.Contains("Syst")){
    if     (Option.Contains("Up"))   SystDir= 1.;
    else if(Option.Contains("Down")) SystDir=-1.;
  }

  bool SiglMuTrig=false, SiglElTrig=false;
  float MinPt=-1, MaxPt=-1, MaxfEta=-1;
  TH2F* HistEff=NULL;
  if(DataYear==2016 && SFKey.Contains("IsoORTkIsoMu24_POGTight")){
    SiglMuTrig=true, MinPt=26., MaxPt=500., MaxfEta=2.4; 
    HistEff = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_IsoMu24_POGTight"];
  }
  else if(DataYear==2017 && SFKey.Contains("IsoMu27_POGTight")){
    SiglMuTrig=true, MinPt=29., MaxPt=1200., MaxfEta=2.4; 
    HistEff = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_IsoMu27_POGTight"];
  }
  else if(DataYear==2018 && SFKey.Contains("IsoMu24_POGTight")){
    SiglMuTrig=true, MinPt=26., MaxPt=1200., MaxfEta=2.4; 
    HistEff = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_IsoMu24_POGTight"];
  }
  else if(DataYear==2016 && SFKey.Contains("Ele27WPTight_POGMVAIsoWP90")){
    SiglElTrig=true, MinPt=35., MaxPt=500., MaxfEta=2.5;
    HistEff = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_Ele27WPTight_POGMVAIsoWP90"];
  }
  else if(DataYear==2017 && SFKey.Contains("Ele32WPTight1OR2_POGMVAIsoWP90")){
    SiglElTrig=true, MinPt=35., MaxPt=500., MaxfEta=2.5;
    HistEff = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_Ele32WPTight1OR2_POGMVAIsoWP90"];
  }
  else if(DataYear==2018 && SFKey.Contains("Ele32WPTight_POGMVAIsoWP90")){
    SiglElTrig=true, MinPt=35., MaxPt=500., MaxfEta=2.5;
    HistEff = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_Ele32WPTight_POGMVAIsoWP90"];
  }

  if(!HistEff){cerr<<"[MCCorrection::MuonTrigger_Eff] No eff file for "<<SFKey<<endl; exit(EXIT_FAILURE);}


  float TriggerEff=0.;
  if(SiglMuTrig){
    for(unsigned int it_m=0; it_m<MuColl.size(); it_m++){
      float pt   = MuColl.at(it_m).Pt();
      float feta = fabs(MuColl.at(it_m).Eta());
      if     (pt<MinPt) return 1.;
      else if(pt>MaxPt) return 1.;
      if     (feta>MaxfEta) return 1.;

      int BinIdx = HistEff->FindBin(pt, feta);
      TriggerEff = HistEff->GetBinContent(BinIdx);
      if(SystDir!=0){ TriggerEff += float(SystDir)*HistEff->GetBinError(BinIdx); }
    }
  }
  else if(SiglElTrig){
    for(unsigned int it_e=0; it_e<EleColl.size(); it_e++){
      float pt   = EleColl.at(it_e).Pt();
      float feta = fabs(EleColl.at(it_e).Eta());
      if     (pt<MinPt) return 1.;
      else if(pt>MaxPt) return 1.;
      if     (feta>MaxfEta) return 1.;

      int BinIdx = HistEff->FindBin(pt, feta);
      TriggerEff = HistEff->GetBinContent(BinIdx);
      if(SystDir!=0){ TriggerEff += float(SystDir)*HistEff->GetBinError(BinIdx); }
    }
  }

  return TriggerEff;
}
