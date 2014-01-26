// $Id: runHgg2013Final_7TeV.C,v 1.6 2013/12/11 19:09:45 bendavid Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TProfile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"

#endif

//--------------------------------------------------------------------------------------------------
void runHgg2013Final_7TeV(const char *fileset    = "0000",
			  const char *skim       = "noskim",
			  //const char *dataset = "s11-h120gg-gf-lv3",
			  //const char *dataset = "r11a-pho-j21-v1",
			  const char *dataset = "s11-h120gg-vh-lv3",
                          //const char *dataset = "r11b-pho-j21-v1",
			  const char *book       = "t2mit/filefi/031",
			  const char *catalogDir = "/home/cmsprod/catalog",
			  const char *outputName = "hgg",
			  int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;
  
  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    sprintf(json, "%s", "~");
    //printf(" JSON file was not properly defined. EXIT!\n");
    //return;
  } 
  
  TString jsonFile = TString("/home/mingyang/cms/json/") + TString(json);
  //TString jsonFile = TString("/home/mingyang/cms/json/") + TString("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt");
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) );
  
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
     sprintf(overlap,"%s", "-1.0");
    //printf(" OVERLAP file was not properly defined. EXIT!\n");
    //return;
  } 

  printf("\n Initialization worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  MCProcessSelectionMod *mcselmod = new MCProcessSelectionMod;
  
  MVASystematicsMod *sysMod = new MVASystematicsMod;
  sysMod->SetMCR9Scale(1.0035, 1.0035);  
  sysMod->SetIsData(isData);
  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/mingyang/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModM = new HLTMod("HLTModM");
  hltModM->AddTrigger("HLT_Mu9");
  hltModM->AddTrigger("HLT_Mu11");
  hltModM->AddTrigger("HLT_Mu15_v1");
  hltModM->SetTrigObjsName("MyHltMuonObjs");
  hltModM->SetAbortIfNotAccepted(kFALSE);

  HLTMod *hltModE = new HLTMod("HLTModE");
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,164237);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",165085,165888);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",165900,166967);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",166968,170053);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170054,170759);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170760,173198);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",173199,178380);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",178381,179889);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,189999);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",190000,999999);
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(isData);

  HLTMod *hltModES = new HLTMod("HLTModES");
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",160000,189999);
  hltModES->AddTrigger("HLT_Ele25_WP80_PFMT40_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_PFMT50_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_v*",190000,999999);
  hltModES->SetTrigObjsName("MyHltElecSObjs");
  hltModES->SetAbortIfNotAccepted(isData);

   HLTMod *hltModP = new HLTMod("HLTModP");
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",160000,161176);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v*",161216,165633);

  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",165970,173198);
  
  hltModP->AddTrigger("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v*",165970,166967);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",167039,173198);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",173236,178380);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",178420,189999);
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);

  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  //  SepPUMod->SetUseAllVerteces(kFALSE);
  // SepPUMod->SetVertexName("OutVtxCiC");
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kTRUE);
   
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (4.0);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetIsMC(!isData);

  GoodPVFilterMod *goodPVFilterModE = new GoodPVFilterMod("GoodPVFilterModE");
  goodPVFilterModE->SetOutputName("GoodVertexesE");
  goodPVFilterModE->SetMinVertexNTracks(0);
  goodPVFilterModE->SetMinNDof         (4.0);
  goodPVFilterModE->SetMaxAbsZ         (24.0);
  goodPVFilterModE->SetMaxRho          (2.0);  
  goodPVFilterModE->SetIsMC(!isData);
  
  GoodPVFilterMod *goodPVFilterModES = new GoodPVFilterMod("GoodPVFilterModES");
  goodPVFilterModES->SetOutputName("GoodVertexesES");
  goodPVFilterModES->SetMinVertexNTracks(0);
  goodPVFilterModES->SetMinNDof         (4.0);
  goodPVFilterModES->SetMaxAbsZ         (24.0);
  goodPVFilterModES->SetMaxRho          (2.0);    
  

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  //-----------------------------------
  // Lepton Selection 
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod -> SetPtMin(20);  
  eleIdMod -> SetEtaMax(2.5);
  eleIdMod -> SetApplyEcalFiducial(true);
  eleIdMod -> SetIDType("Hgg_LeptonTag_2012IdHCP");  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml")))); 
  eleIdMod -> SetWhichVertex(-1);
  eleIdMod -> SetD0Cut(0.02);
  eleIdMod -> SetDZCut(0.2); //h  
  eleIdMod -> SetIsoType("PFIso_HggLeptonTag2012HCP"); //h
  eleIdMod -> SetOutputName("HggLeptonTagElectrons");
  eleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  eleIdMod -> SetPFNoPileUpName("pfnopileupcands");
  eleIdMod -> SetInvertNExpectedHitsInnerCut(kFALSE);
  eleIdMod -> SetNExpectedHitsInnerCut(1);   
  eleIdMod -> SetApplyConversionFilterType1(kTRUE);
  eleIdMod -> SetPVName(Names::gkPVBeamSpotBrn);   

  ElectronIDMod* softEleIdMod = new ElectronIDMod;
  softEleIdMod -> SetPtMin(10);  
  softEleIdMod -> SetEtaMax(2.5);
  softEleIdMod -> SetApplyEcalFiducial(true);
  softEleIdMod -> SetIDType("Hgg_LeptonTag_2012IdHCP");  
  softEleIdMod -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  softEleIdMod -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  softEleIdMod -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  softEleIdMod -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  softEleIdMod -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  softEleIdMod -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml")))); 
  softEleIdMod -> SetWhichVertex(-1);
  softEleIdMod -> SetD0Cut(0.02);
  softEleIdMod -> SetDZCut(0.2); //h  
  softEleIdMod -> SetIsoType("PFIso_HggLeptonTag2012HCP"); //h
  softEleIdMod -> SetOutputName("HggLeptonTagSoftElectrons");
  softEleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  softEleIdMod -> SetPFNoPileUpName("pfnopileupcands");
  softEleIdMod -> SetInvertNExpectedHitsInnerCut(kFALSE);
  softEleIdMod -> SetNExpectedHitsInnerCut(1);   
  softEleIdMod -> SetApplyConversionFilterType1(kTRUE);
  softEleIdMod -> SetPVName(Names::gkPVBeamSpotBrn);
   
  MuonIDMod* muonIdMod = new MuonIDMod;
  // base kinematics
  muonIdMod -> SetPtMin(20.);
  muonIdMod -> SetEtaCut(2.4);
  // base ID
  //muonIdMod -> SetIDType("Tight");
  muonIdMod -> SetIDType("muonPOG2012CutBasedIDTight");
  muonIdMod -> SetWhichVertex(-1); // this is a 'hack'.. but hopefully good enough...
  muonIdMod -> SetD0Cut(0.2);
  muonIdMod -> SetDZCut(0.5);
  muonIdMod -> SetIsoType("PFIsoBetaPUCorrected"); //h
  muonIdMod -> SetPFIsoCut(0.2); //h
  muonIdMod -> SetOutputName("HggLeptonTagMuons");
  muonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  muonIdMod -> SetPFPileUpName("pfpileupcands");
  muonIdMod -> SetPVName(Names::gkPVBeamSpotBrn); 

 MuonIDMod* softMuonIdMod = new MuonIDMod;
  // base kinematics
  softMuonIdMod -> SetPtMin(10.);
  softMuonIdMod -> SetEtaCut(2.4);
  // base ID
  //softMuonIdMod -> SetIDType("Tight");
  softMuonIdMod -> SetIDType("muonPOG2012CutBasedIDTight");
  softMuonIdMod -> SetWhichVertex(-1); // this is a 'hack'.. but hopefully good enough...
  softMuonIdMod -> SetD0Cut(0.2);
  softMuonIdMod -> SetDZCut(0.5);
  softMuonIdMod -> SetIsoType("PFIsoBetaPUCorrected"); //h
  softMuonIdMod -> SetPFIsoCut(0.2); //h
  softMuonIdMod -> SetOutputName("HggLeptonTagSoftMuons");
  softMuonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  softMuonIdMod -> SetPFPileUpName("pfpileupcands");
  softMuonIdMod -> SetPVName(Names::gkPVBeamSpotBrn); 
  
  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  
  
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
   if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/2011ReRecoJEC_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data())); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    
  
  Bool_t excludedoubleprompt = kFALSE;
  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludedoubleprompt = kTRUE;
  }
  
  if (TString(dataset).Contains("-qcd")) {
    excludedoubleprompt = kTRUE;
  }
  
  Bool_t is25 = kFALSE;
  if (TString(book).Contains("025")) is25 = kTRUE;
  
  TString encorrfilename = std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/21Jun2012_7TeV-step2-invMass_SC_regrCorrSemiPar7TeVtrainV8_pho-loose-Et_25-noPF-HggRunEtaR9.dat")).Data());
  
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetRegressionVersion(8);
  photreg->SetRegressionWeights(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/regweights_v8_7TeV_forest_ph.root")).Data()));
  photreg->SetOutputName("GoodPhotonsRegr");
  //photreg->SetApplyShowerRescaling(kTRUE);
  photreg->SetIsData(isData);

  
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetR9Rescale(true, 0.000854, 1.00153, 0.001231,1.0005);
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetOutputVtxName("OutVtxCiC");        
  photcic->SetPhotonSelType("CiCPFSelection");
  photcic->SetVertexSelType("CiCMVASelection");
  photcic->SetUseSingleLegConversions(kFALSE);
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kTRUE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcic->SetMCSmearFactors(0.0068, 0.0096, 0.0101, 0.0185, 0.0158, 0.0185, 0.0201, 0.0183);
  
  photcic->SetMCStochasticPivot(6.60119,6.7276,6.52397,6.73261);
  photcic->SetMCStochasticRho(0.0068, 0.0096, 0.0101, 0.0185);
  photcic->SetMCStochasticPhi(TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.);
  photcic->SetStochasticSmear(false);  
  
  photcic->AddEnCorrFromFile(encorrfilename);

  //-----------------------------------------------------------------------------------------------------------------
  
  //  photcic->SetDoMCR9Scaling(kTRUE);
  //  photcic->SetMCR9Scale(1.0035, 1.0035);
  photcic->SetDoShowerShapeScaling(kFALSE); 
  photcic->SetShowerShapeType("2011ShowerShape");
  //photcic->SetDoMCErrScaling(kTRUE);
  //photcic->SetMCErrScale(1.07, 1.045); 
  //photcic->SetMCErrScale(1, 1); //ming:scale(sigE/E)
  photcic->SetJetsName(jetCorr->GetOutputName());    
  //photcic->SetRescaledBeamspotWidth(5.0);
  photcic->SetIsData(isData);
  photcic->SetApplyLeptonTag(kTRUE);
  photcic->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
  photcic->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");  
  photcic->SetIdMVAType("2013FinalIdMVA_7TeV");
  
  PhotonPairSelector         *photcicnoeleveto = new PhotonPairSelector("PhotonPairSelectorCiCInvertEleVeto");
  photcicnoeleveto->SetR9Rescale(true, 0.000854, 1.00153, 0.001231,1.0005);
  photcicnoeleveto->SetOutputName("GoodPhotonsCICNoEleVeto");
  photcicnoeleveto->SetOutputVtxName("OutVtxCiCInvertEleVeto");      
  photcicnoeleveto->SetPhotonSelType("CiCPFSelection");
  photcicnoeleveto->SetVertexSelType("CiCMVASelection");
  photcicnoeleveto->SetUseSingleLegConversions(kFALSE);
  photcicnoeleveto->DoMCSmear(kTRUE);
  photcicnoeleveto->DoDataEneCorr(kTRUE);
  photcicnoeleveto->SetPhotonsFromBranch(kFALSE);
  photcicnoeleveto->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcicnoeleveto->SetMCSmearFactors(0.0068, 0.0096, 0.0101, 0.0185, 0.0158, 0.0185, 0.0201, 0.0183);

  photcicnoeleveto->SetMCStochasticPivot(6.60119,6.7276,6.52397,6.73261);
  photcicnoeleveto->SetMCStochasticRho(0.0068, 0.0096, 0.0101, 0.0185);
  photcicnoeleveto->SetMCStochasticPhi(TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.);
  photcicnoeleveto->SetStochasticSmear(false);    
  
  photcicnoeleveto->AddEnCorrFromFile(encorrfilename);


  //-----------------------------------------------------------------------------------------------------------------
  
  //photcicnoeleveto->SetDoMCR9Scaling(kTRUE);
  //photcicnoeleveto->SetMCR9Scale(1.0035, 1.0035);
  photcicnoeleveto->SetDoShowerShapeScaling(kFALSE);
  photcicnoeleveto->SetShowerShapeType("2011ShowerShape");
  //photcicnoeleveto->SetDoMCErrScaling(kTRUE);
  //photcicnoeleveto->SetMCErrScale(1.07, 1.045);    
  //photcicnoeleveto->SetMCErrScale(1, 1);    
  photcicnoeleveto->SetApplyEleVeto(kFALSE);
  photcicnoeleveto->SetInvertElectronVeto(kTRUE);
  photcicnoeleveto->SetJetsName(jetCorr->GetOutputName());  
  //photcicnoeleveto->SetRescaledBeamspotWidth(5.0);
  photcicnoeleveto->SetIsData(isData);
  photcicnoeleveto->SetApplyLeptonTag(kTRUE);
  photcicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
  photcicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");  
  photcicnoeleveto->SetIdMVAType("2013FinalIdMVA_7TeV");
 
  PhotonPairSelector         *photpresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photpresel->SetOutputName("GoodPhotonsPresel");
  photpresel->SetPhotonSelType("MITPFSelection");//diff 2012
  //photpresel->SetPhotonSelType("MITSelection");
  photpresel->SetVertexSelType("CiCMVASelection");//diff 2012
  photpresel->SetUseSingleLegConversions(kFALSE);
  photpresel->SetIdMVAType("2013FinalIdMVA_7TeV");//diff 2012
  photpresel->DoMCSmear(kTRUE);
  photpresel->DoDataEneCorr(kTRUE);
  photpresel->SetPhotonsFromBranch(kFALSE);
  photpresel->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpresel->SetMCSmearFactors(0.0068, 0.0096, 0.0101, 0.0185, 0.0158, 0.0185, 0.0201, 0.0183);

  photpresel->SetMCStochasticPivot(6.60119,6.7276,6.52397,6.73261);
  photpresel->SetMCStochasticRho(0.0068, 0.0096, 0.0101, 0.0185);
  photpresel->SetMCStochasticPhi(TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.);
  photpresel->SetStochasticSmear(false);      
  
  photpresel->AddEnCorrFromFile(encorrfilename);

  
  //-----------------------------------------------------------------------------------------------------------------
  //photpresel->SetDoMCR9Scaling(kTRUE);
  photpresel->SetDoShowerShapeScaling(kFALSE);
  photpresel->SetShowerShapeType("2011ShowerShape");
  //photpresel->SetDoMCErrScaling(kTRUE);
  //photpresel->SetMCErrScale(1.07, 1.045); 
  //photpresel->SetMCErrScale(1, 1); 
  photpresel->SetJetsName(jetCorr->GetOutputName()); 
  //photpresel->SetRescaledBeamspotWidth(5.0);
  photpresel->SetIsData(isData);
  photpresel->SetApplyLeptonTag(kTRUE);
  photpresel->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
  photpresel->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");  
 
  PhotonPairSelector         *photpreselinverteleveto = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputName("GoodPhotonsPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputVtxName("OutVtxPreselInvertEleVeto");    
  photpreselinverteleveto->SetPhotonSelType("MITPFSelection");
  photpreselinverteleveto->SetIdMVAType("2013FinalIdMVA_7TeV");
  //photpreselinverteleveto->SetVertexSelType("CiCMVA2012Selection");//ming change
  photpreselinverteleveto->SetUseSingleLegConversions(kFALSE);
  photpreselinverteleveto->DoMCSmear(kTRUE);
  photpreselinverteleveto->DoDataEneCorr(kTRUE);
  photpreselinverteleveto->SetPhotonsFromBranch(kFALSE);
  photpreselinverteleveto->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpreselinverteleveto->SetMCSmearFactors(0.0068, 0.0096, 0.0101, 0.0185, 0.0158, 0.0185, 0.0201, 0.0183);
  
  photpreselinverteleveto->SetMCStochasticPivot(6.60119,6.7276,6.52397,6.73261);
  photpreselinverteleveto->SetMCStochasticRho(0.0068, 0.0096, 0.0101, 0.0185);
  photpreselinverteleveto->SetMCStochasticPhi(TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.,TMath::Pi()/2.);
  photpreselinverteleveto->SetStochasticSmear(false);     
  
  photpreselinverteleveto->AddEnCorrFromFile(encorrfilename);
  
  
  //-----------------------------------------------------------------------------------------------------------------
  photpreselinverteleveto->SetShowerShapeType("2011ShowerShape");
  photpreselinverteleveto->SetDoShowerShapeScaling(kFALSE);  
  //photpreselinverteleveto->SetDoMCErrScaling(kTRUE);
  //photpreselinverteleveto->SetMCErrScale(1.07, 1.045);  
  //photpreselinverteleveto->SetMCErrScale(1, 1);  
  photpreselinverteleveto->SetApplyEleVeto(kFALSE);
  photpreselinverteleveto->SetInvertElectronVeto(kTRUE);
  photpreselinverteleveto->SetJetsName(jetCorr->GetOutputName());  
  //photpreselinverteleveto->SetRescaledBeamspotWidth(5.0);   
  photpreselinverteleveto->SetIsData(isData);
  photpreselinverteleveto->SetApplyLeptonTag(kTRUE);
  photpreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
  photpreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");    
  photpreselinverteleveto->SetLeadingPtMin(20.);
  photpreselinverteleveto->SetTrailingPtMin(20.);

   PhotonPairSelector         *photpreselnosmear = new PhotonPairSelector("PhotonPairSelectorPreselNoSmear");
   photpreselnosmear->SetOutputName("GoodPhotonsPreselNoSmear");
   photpreselnosmear->SetPhotonSelType("MITPFSelection");
   photpreselnosmear->SetVertexSelType("CiCMVASelection");
   photpreselnosmear->SetUseSingleLegConversions(kFALSE);  
   photpreselnosmear->SetIdMVAType("2013FinalIdMVA_7TeV");
   photpreselnosmear->SetShowerShapeType("2011ShowerShape");
   photpreselnosmear->SetDoShowerShapeScaling(kFALSE);  
   photpreselnosmear->SetPhotonsFromBranch(kFALSE);
   photpreselnosmear->SetInputPhotonsName(photreg->GetOutputName());
   photpreselnosmear->SetMCSmearFactors(0.0068, 0.0096, 0.0101, 0.0185, 0.0158, 0.0185, 0.0201, 0.0183);
   photpreselnosmear->SetJetsName(jetCorr->GetOutputName());  
   photpreselnosmear->SetOutputVtxName("OutVtxNoSmear");  
   photpreselnosmear->SetLeadingPtMin(30.);
   photpreselnosmear->SetTrailingPtMin(22.);
   //photpreselnosmear->SetRescaledBeamspotWidth(5.0);      
   photpreselnosmear->SetIsData(isData);  
   photpreselnosmear->SetApplyLeptonTag(kTRUE);
   photpreselnosmear->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   photpreselnosmear->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   photpreselnosmear->DoMCSmear(kFALSE);
   photpreselnosmear->DoMCEneSmear(kFALSE);
   photpreselnosmear->DoEneErrSmear(kTRUE);
   photpreselnosmear->DoDataEneCorr(kFALSE);

   PhotonPairSelector         *photcicnosmear = new PhotonPairSelector("PhotonPairSelectorCiCNoSmear");
   photcicnosmear->SetR9Rescale(true, 0.000854, 1.00153, 0.001231,1.0005);
   photcicnosmear->SetOutputName("GoodPhotonsCICNoSmear");
   photcicnosmear->SetOutputVtxName("OutVtxCiCNoSmear");        
   photcicnosmear->SetPhotonSelType("CiCPFSelection");
   photcicnosmear->SetVertexSelType("CiCMVASelection");
   photcicnosmear->SetUseSingleLegConversions(kFALSE);
   photcicnosmear->DoMCSmear(kFALSE);
   photcicnosmear->DoDataEneCorr(kFALSE);
   photcicnosmear->SetPhotonsFromBranch(kFALSE);
   photcicnosmear->SetInputPhotonsName(photreg->GetOutputName());
   //  photcicnosmear->SetDoMCR9Scaling(kTRUE);
   //  photcicnosmear->SetMCR9Scale(1.0035, 1.0035);
   photcicnosmear->SetDoShowerShapeScaling(kFALSE); 
   photcicnosmear->SetShowerShapeType("2011ShowerShape");
   //photcicnosmear->SetDoMCErrScaling(kTRUE);
   //photcicnosmear->SetMCErrScale(1.07, 1.045); 
   //photcicnosmear->SetMCErrScale(1, 1); //ming:scale(sigE/E)
   photcicnosmear->SetJetsName(jetCorr->GetOutputName());    
   //photcicnosmear->SetRescaledBeamspotWidth(5.0);
   photcicnosmear->SetIsData(isData);
   photcicnosmear->SetApplyLeptonTag(kTRUE);
   photcicnosmear->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   photcicnosmear->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");  
   
   PhotonPairSelector         *photpreselinvertelevetonosmear = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVetoNoSmear");
   photpreselinvertelevetonosmear->SetOutputName("GoodPhotonsPreselInvertEleVetoNoSmear");
   photpreselinvertelevetonosmear->SetPhotonSelType("MITPFSelectionNoEcal");
   //photpreselinvertelevetonosmear->SetVertexSelType("CiCMVA2012Selection");//ming change
   photpreselinvertelevetonosmear->SetUseSingleLegConversions(kFALSE);  
   photpreselinvertelevetonosmear->SetIdMVAType("2013FinalIdMVA_7TeV");
   photpreselinvertelevetonosmear->SetShowerShapeType("2011ShowerShape");
   photpreselinvertelevetonosmear->SetDoShowerShapeScaling(kFALSE);  
   photpreselinvertelevetonosmear->SetPhotonsFromBranch(kFALSE);
   photpreselinvertelevetonosmear->SetInputPhotonsName(photreg->GetOutputName());
   photpreselinvertelevetonosmear->SetJetsName(jetCorr->GetOutputName());  
   photpreselinvertelevetonosmear->SetOutputVtxName("OutVtxPreselInvertEleVetoNoSmear");  
   //photpreselinvertelevetonosmear->SetLeadingPtMin(30.);
   photpreselinvertelevetonosmear->SetLeadingPtMin(20.);//ming
   //photpreselinvertelevetonosmear->SetTrailingPtMin(22.);
   photpreselinvertelevetonosmear->SetTrailingPtMin(20.);//ming
   //photpreselinvertelevetonosmear->SetRescaledBeamspotWidth(5.0);      
   photpreselinvertelevetonosmear->SetIsData(isData);  
   photpreselinvertelevetonosmear->SetApplyLeptonTag(kTRUE);
   photpreselinvertelevetonosmear->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   photpreselinvertelevetonosmear->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   photpreselinvertelevetonosmear->DoMCSmear(kFALSE);
   photpreselinvertelevetonosmear->DoDataEneCorr(kFALSE);
   photpreselinvertelevetonosmear->SetApplyEleVeto(kFALSE);
   photpreselinvertelevetonosmear->SetInvertElectronVeto(kTRUE);
   
   PhotonTreeWriter *phottreecic = new PhotonTreeWriter("PhotonTreeWriterCiC");
   phottreecic->SetR9Rescale(true, 0.000854, 1.00153, 0.001231,1.0005);
   phottreecic->SetPhotonsFromBranch(kFALSE);
   phottreecic->SetInputPhotonsName(photcic->GetOutputName());
   phottreecic->SetEnableJets(kTRUE);
   phottreecic->SetApplyJetId(kFALSE);
   phottreecic->SetPFJetsFromBranch(kFALSE);
   phottreecic->SetPFJetName(jetCorr->GetOutputName());
   phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecic->SetIsData(isData);
   phottreecic->SetIsCutBased(kTRUE);
   if (is25) phottreecic->SetEnablePFPhotons(kFALSE);
   phottreecic->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreecic->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreecic->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreecic->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreecic->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreecic->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   //phottreecic->SetApplyLeptonTag(kTRUE);
   //phottreecic->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   //phottreecic->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   //phottreecic->SetApplyVBFTag(kTRUE);
   phottreecic->SetApplyPFMetCorr(kTRUE);
   phottreecic->SetApplyVHLepTag(kTRUE);
   phottreecic->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreecic->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreecic->SetApplyVHHadTag(kTRUE);
   phottreecic->SetApplyVBFTag(kTRUE);
   phottreecic->SetApplyTTHTag(kTRUE);
   
   PhotonTreeWriter *phottreecicnoeleveto = new PhotonTreeWriter("PhotonTreeWriterCiCInvertEleVeto");
   phottreecicnoeleveto->SetR9Rescale(true, 0.000854, 1.00153, 0.001231,1.0005);
   phottreecicnoeleveto->SetPhotonsFromBranch(kFALSE);
   phottreecicnoeleveto->SetInputPhotonsName(photcicnoeleveto->GetOutputName());
   phottreecicnoeleveto->SetEnableJets(kTRUE);
   phottreecicnoeleveto->SetApplyJetId(kFALSE);  
   phottreecicnoeleveto->SetPFJetsFromBranch(kFALSE);
   phottreecicnoeleveto->SetPFJetName(jetCorr->GetOutputName());
   phottreecicnoeleveto->SetApplyElectronVeto(kFALSE);
   phottreecicnoeleveto->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicnoeleveto->SetIsData(isData);  
   phottreecicnoeleveto->SetIsCutBased(kTRUE);
   if (is25) phottreecicnoeleveto->SetEnablePFPhotons(kFALSE);
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   //phottreecicnoeleveto->SetApplyLeptonTag(kTRUE);
   //phottreecicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   //phottreecicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   //phottreecicnoeleveto->SetApplyVBFTag(kTRUE);
   phottreecicnoeleveto->SetApplyPFMetCorr(kTRUE);
   phottreecicnoeleveto->SetApplyVHLepTag(kTRUE);
   phottreecicnoeleveto->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreecicnoeleveto->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreecicnoeleveto->SetApplyVHHadTag(kTRUE);
   phottreecicnoeleveto->SetApplyVBFTag(kTRUE);
   phottreecicnoeleveto->SetApplyTTHTag(kTRUE);
   
   PhotonTreeWriter *phottreepresel = new PhotonTreeWriter("PhotonTreeWriterPresel");
   phottreepresel->SetPhotonsFromBranch(kFALSE);
   phottreepresel->SetInputPhotonsName(photpresel->GetOutputName());
   phottreepresel->SetEnableJets(kTRUE);
   phottreepresel->SetApplyJetId(kFALSE);    
   phottreepresel->SetPFJetsFromBranch(kFALSE);
   phottreepresel->SetPFJetName(jetCorr->GetOutputName());  
   phottreepresel->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepresel->SetIsData(isData);  
   phottreepresel->SetIsCutBased(kFALSE);
   if (is25) phottreepresel->SetEnablePFPhotons(kFALSE);
   phottreepresel->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreepresel->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreepresel->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreepresel->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreepresel->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreepresel->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   phottreepresel->SetDoSynching(kTRUE);
   // phottreepresel->SetApplyLeptonTag(kTRUE);
   // phottreepresel->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   // phottreepresel->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   // phottreepresel->SetApplyVBFTag(kTRUE);
    phottreepresel->SetApplyPFMetCorr(kTRUE);
    phottreepresel->SetApplyVHLepTag(kTRUE);
    phottreepresel->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
    phottreepresel->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
    phottreepresel->SetApplyVHHadTag(kTRUE);
    phottreepresel->SetApplyVBFTag(kTRUE);
    phottreepresel->SetApplyTTHTag(kTRUE);
    
   PhotonTreeWriter *phottreepreselinverteleveto = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVeto");
   phottreepreselinverteleveto->SetPhotonsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetInputPhotonsName(photpreselinverteleveto->GetOutputName());
   phottreepreselinverteleveto->SetEnableJets(kTRUE);
   phottreepreselinverteleveto->SetApplyJetId(kFALSE);    
   phottreepreselinverteleveto->SetPFJetsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselinverteleveto->SetApplyElectronVeto(kFALSE);  
   phottreepreselinverteleveto->SetExcludeDoublePrompt(excludedoubleprompt);    
   phottreepreselinverteleveto->SetIsData(isData);
   phottreepreselinverteleveto->SetIsCutBased(kFALSE); 
   if (is25) phottreepreselinverteleveto->SetEnablePFPhotons(kFALSE);  
   phottreepreselinverteleveto->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreepreselinverteleveto->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreepreselinverteleveto->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreepreselinverteleveto->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreepreselinverteleveto->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreepreselinverteleveto->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   //  phottreepreselinverteleveto->SetApplyLeptonTag(kTRUE);
   //  phottreepreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   //  phottreepreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   //  phottreepreselinverteleveto->SetApplyVBFTag(kTRUE);
   phottreepreselinverteleveto->SetApplyPFMetCorr(kTRUE);
   phottreepreselinverteleveto->SetApplyVHLepTag(kTRUE);
   phottreepreselinverteleveto->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreepreselinverteleveto->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreepreselinverteleveto->SetApplyVHHadTag(kTRUE);
   phottreepreselinverteleveto->SetApplyVBFTag(kTRUE);
   phottreepreselinverteleveto->SetApplyTTHTag(kTRUE);
  
   PhotonTreeWriter *phottreepreselnosmear = new PhotonTreeWriter("PhotonTreeWriterPreselNoSmear");
   phottreepreselnosmear->SetPhotonsFromBranch(kFALSE);
   phottreepreselnosmear->SetInputPhotonsName(photpreselnosmear->GetOutputName());
   phottreepreselnosmear->SetEnableJets(kTRUE);
   phottreepreselnosmear->SetApplyJetId(kFALSE);  
   phottreepreselnosmear->SetPFJetsFromBranch(kFALSE);
   phottreepreselnosmear->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselnosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepreselnosmear->SetIsData(isData);    
   phottreepreselnosmear->SetIsCutBased(kFALSE);
   if (is25) phottreepreselnosmear->SetEnablePFPhotons(kFALSE); 
   phottreepreselnosmear->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreepreselnosmear->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreepreselnosmear->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreepreselnosmear->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreepreselnosmear->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreepreselnosmear->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   phottreepreselnosmear->SetDoSynching(kTRUE);
   //  phottreepreselnosmear->SetApplyLeptonTag(kTRUE);
   //  phottreepreselnosmear->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   //  phottreepreselnosmear->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   //  phottreepreselnosmear->SetApplyVBFTag(kTRUE);
   phottreepreselnosmear->SetApplyPFMetCorr(kTRUE);
   phottreepreselnosmear->SetApplyVHLepTag(kTRUE);
   phottreepreselnosmear->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreepreselnosmear->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreepreselnosmear->SetApplyVHHadTag(kTRUE);
   phottreepreselnosmear->SetApplyVBFTag(kTRUE);
   phottreepreselnosmear->SetApplyTTHTag(kTRUE);

   PhotonTreeWriter *phottreecicnosmear = new PhotonTreeWriter("PhotonTreeWriterCiCNoSmear");
   phottreecicnosmear->SetR9Rescale(true, 0.000854, 1.00153, 0.001231,1.0005);
   phottreecicnosmear->SetPhotonsFromBranch(kFALSE);
   phottreecicnosmear->SetInputPhotonsName(photcicnosmear->GetOutputName());
   phottreecicnosmear->SetEnableJets(kTRUE);
   phottreecicnosmear->SetApplyJetId(kFALSE);
   phottreecicnosmear->SetPFJetsFromBranch(kFALSE);
   phottreecicnosmear->SetPFJetName(jetCorr->GetOutputName());
   phottreecicnosmear->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicnosmear->SetIsData(isData);
   phottreecicnosmear->SetIsCutBased(kTRUE);
   if (is25) phottreecicnosmear->SetEnablePFPhotons(kFALSE);
   phottreecicnosmear->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreecicnosmear->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreecicnosmear->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreecicnosmear->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreecicnosmear->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreecicnosmear->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   //  phottreecicnosmear->SetApplyLeptonTag(kTRUE);
   //  phottreecicnosmear->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   //  phottreecicnosmear->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   //  phottreecicnosmear->SetApplyVBFTag(kTRUE);
   phottreecicnosmear->SetApplyPFMetCorr(kTRUE);
   phottreecicnosmear->SetApplyVHLepTag(kTRUE);
   phottreecicnosmear->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreecicnosmear->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreecicnosmear->SetApplyVHHadTag(kTRUE);
   phottreecicnosmear->SetApplyVBFTag(kTRUE);
   phottreecicnosmear->SetApplyTTHTag(kTRUE);

   PhotonTreeWriter *phottreepreselinvertelevetonosmear = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVetoNoSmear");
   phottreepreselinvertelevetonosmear->SetPhotonsFromBranch(kFALSE);
   phottreepreselinvertelevetonosmear->SetInputPhotonsName(photpreselinvertelevetonosmear->GetOutputName());
   phottreepreselinvertelevetonosmear->SetEnableJets(kTRUE);
   phottreepreselinvertelevetonosmear->SetApplyJetId(kFALSE);  
   phottreepreselinvertelevetonosmear->SetPFJetsFromBranch(kFALSE);
   phottreepreselinvertelevetonosmear->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselinvertelevetonosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepreselinvertelevetonosmear->SetIsData(isData);  
   phottreepreselinvertelevetonosmear->SetIsCutBased(kFALSE);  
   if (is25) phottreepreselinvertelevetonosmear->SetEnablePFPhotons(kFALSE); 
   phottreepreselinvertelevetonosmear->SetApplyElectronVeto(kFALSE); 
  phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   // phottreepreselinvertelevetonosmear->SetApplyLeptonTag(kTRUE);
   // phottreepreselinvertelevetonosmear->SetLeptonTagElectronsName("HggLeptonTagSoftElectrons");
   // phottreepreselinvertelevetonosmear->SetLeptonTagMuonsName("HggLeptonTagSoftMuons");
   // phottreepreselinvertelevetonosmear->SetApplyVBFTag(kTRUE);
  phottreepreselinvertelevetonosmear->SetApplyPFMetCorr(kTRUE);
  phottreepreselinvertelevetonosmear->SetApplyVHLepTag(kTRUE);
  phottreepreselinvertelevetonosmear->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
  phottreepreselinvertelevetonosmear->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
  phottreepreselinvertelevetonosmear->SetApplyVHHadTag(kTRUE);
  phottreepreselinvertelevetonosmear->SetApplyVBFTag(kTRUE);
  phottreepreselinvertelevetonosmear->SetApplyTTHTag(kTRUE);
   
   PhotonIDMod         *photidpresel = new PhotonIDMod("PhotonIDModPresel");
   photidpresel->SetPtMin(25.0);
   photidpresel->SetOutputName("GoodPhotonsPreselid");
   photidpresel->SetIDType("MITPFSelection");
   photidpresel->SetApplyElectronVeto(kTRUE);
   photidpresel->SetIsData(isData);
   
  PhotonIDMod         *photidpreselinvert = new PhotonIDMod("PhotonIDModPreselInvert");
  photidpreselinvert->SetPtMin(25.0);
  photidpreselinvert->SetOutputName("GoodPhotonsPreselidInvert");
  photidpreselinvert->SetIDType("MITPFSelection");
  photidpreselinvert->SetApplyElectronVeto(kFALSE);
  photidpreselinvert->SetInvertElectronVeto(kTRUE);
  photidpreselinvert->SetIsData(isData);  
  
  PhotonTreeWriter *phottreesingle = new PhotonTreeWriter("PhotonTreeWriterSingle");
  phottreesingle->SetWriteDiphotonTree(kFALSE);
  phottreesingle->SetPhotonsFromBranch(kFALSE);
  phottreesingle->SetInputPhotonsName(photidpresel->GetOutputName());  
  phottreesingle->SetEnableJets(kTRUE);
  phottreesingle->SetPFJetsFromBranch(kFALSE);
  phottreesingle->SetPFJetName(jetCorr->GetOutputName());      
  //phottreesingle->SetBeamspotWidth(5.0);
  phottreesingle->SetIsData(isData);
  if (is25) phottreesingle->SetEnablePFPhotons(kFALSE);  
  
  PhotonTreeWriter *phottreesingleinvert = new PhotonTreeWriter("PhotonTreeWriterSingleInvert");
  phottreesingleinvert->SetWriteDiphotonTree(kTRUE);
  phottreesingleinvert->SetPhotonsFromBranch(kFALSE);
  phottreesingleinvert->SetInputPhotonsName(photidpreselinvert->GetOutputName());  
  phottreesingleinvert->SetEnableJets(kTRUE);
  phottreesingleinvert->SetPFJetsFromBranch(kFALSE);
  phottreesingleinvert->SetPFJetName(pubJetOpen->GetOutputName());
  phottreesingleinvert->SetApplyElectronVeto(kFALSE);
  //phottreesingleinvert->SetBeamspotWidth(5.0);
  phottreesingleinvert->SetIsData(isData);  
  if (is25) phottreesingleinvert->SetEnablePFPhotons(kFALSE);    
  

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
  
  // high level trigger is always first
  //mcselmod         ->Add(hltModES);
  //mcselmod         ->Add(hltModE);
  
  if (!TString(dataset).Contains("meridiani")) {      
    mcselmod         ->Add(hltModP);
    hltModP         ->Add(goodPVFilterMod);
  }
  else {
    mcselmod->Add(goodPVFilterMod);
  }
  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  //goodPVFilterMod ->Add(muonId);
  goodPVFilterMod->Add(photreg);
  photreg->Add(pubJet);
  pubJet->Add(jetCorr);
  
  // simple object id modules
  //goodPVFilterModE -> Add(pubJetOpen);
  //pubJetOpen       -> Add(photidpreselinvert);
  //goodPVFilterModES ->Add(elecIdS);
 
  jetCorr          ->Add(SepPUMod); 
  SepPUMod         ->Add(eleIdMod);
  eleIdMod         ->Add(softEleIdMod);
  softEleIdMod     ->Add(softMuonIdMod);
  softMuonIdMod    ->Add(muonIdMod);
     
  muonIdMod          ->Add(photcic);
  muonIdMod          ->Add(photcicnoeleveto);  
  muonIdMod          ->Add(photpresel);  
  muonIdMod          ->Add(photpreselinverteleveto);  
  muonIdMod          ->Add(photpreselnosmear);  
  muonIdMod          ->Add(photpreselinvertelevetonosmear);  
  muonIdMod          ->Add(photcicnosmear);  

  photcic         ->Add(phottreecic);
  photcicnoeleveto       ->Add(phottreecicnoeleveto);
  photpresel    ->Add(phottreepresel);
  photpreselinverteleveto    ->Add(phottreepreselinverteleveto);
  photpreselnosmear    ->Add(phottreepreselnosmear);
  photpreselinvertelevetonosmear    ->Add(phottreepreselinvertelevetonosmear);
  photcicnosmear    ->Add(phottreecicnosmear);

  //jetCorr          ->Add(photidpresel);
  //photidpresel       ->Add(phottreesingle);
  
  //photidpreselinvert       ->Add(phottreesingleinvert);  
  //elecIdS->Add(phottreeES);
  

  //TFile::SetOpenTimeout(0);
  //TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
  TFile::SetReadaheadSize(128*1024*1024);
  
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  if (TString(dataset).Contains("meridiani")) {    
    ana->SetUseHLT(kFALSE);
  }

 
  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  Bool_t caching = kTRUE;
  //Bool_t caching = kFALSE;
  //if (TString(dataset).Contains("s11-h")) bookstr.ReplaceAll("local","t2mit");
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset,caching);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset,caching);

  ana->AddDataset(d);
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12a-pho-j13-v1/385B77DE-58D0-E111-B925-001E67396928.root");//79737729
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/FE308E0A-23D2-E111-8B2D-00266CFAE7D0.root");//871378986
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/54E239DC-2725-E211-A52C-003048D373F6.root");//528937923
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/EC80B121-7E30-E211-BFFF-5404A640A63D.root");//261543921
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/66CE9443-472C-E211-B08F-001D09F2437B.root");//483561562
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/F03B81B2-0B29-E211-B28F-BCAEC518FF67.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/A6EC91DD-CB39-E211-8897-BCAEC53296F8.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/64F7EAB0-CBEA-E111-B8AD-003048F117B6.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/0615B16F-ECDA-E111-A266-003048F117B6.root"); 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12a-pho-j13-v1/385B77DE-58D0-E111-B925-001E67396928.root"); 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/068B257D-81D0-E111-A835-5404A63886B4.root");//369441614 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/5443DBCB-562C-E211-8042-00237DDBE0E2.root");//333643114
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/B2698ACF-843A-E211-8D9E-001D09F244DE.root");//89022540
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/1E050D78-6DD4-E111-81D2-00266CFAE8D0.root");//223740859
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/2CD64E47-79D3-E111-90C5-00A0D1EE8E60.root");//541603559
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/147EF338-B2D3-E111-B1A6-00266CF256CC.root");//8983064 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/8843BD9F-DC3E-E211-9E2F-003048678110.root");//876316897 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/CAA12077-F9F6-E111-8DF2-001D09F28D54.root");
 
  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  //ana->SetCacheSize(64*1024*1024);
  ana->SetCacheSize(0);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}
