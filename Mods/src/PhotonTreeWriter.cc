#include <time.h>
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PFMetCorrectionTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "TDataMember.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>

using namespace mithep;

ClassImp(mithep::PhotonTreeWriter)
  templateClassImp(mithep::PhotonTreeWriterPhoton)//ming: what's this?
ClassImp(mithep::PhotonTreeWriterDiphotonEvent)

//--------------------------------------------------------------------------------------------------
PhotonTreeWriter::PhotonTreeWriter(const char *name, const char *title) : 
  // Base Module...
  BaseMod                 (name,title),
  // define all the Branches to load
  fPhotonBranchName       (Names::gkPhotonBrn),
  fPFPhotonName           ("PFPhotons"),
  fElectronName           (Names::gkElectronBrn),
  fGoodElectronName       (Names::gkElectronBrn),  
  fConversionName         (Names::gkMvfConversionBrn),  
  fPFConversionName       ("PFPhotonConversions"),  
  fTrackBranchName        (Names::gkTrackBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fBeamspotName           (Names::gkBeamSpotBrn),
  fPFCandName             (Names::gkPFCandidatesBrn),
  fPFNoPileUpName         ("PFNoPileUp"),
  fPFPileUpName           ("PFPileUp"),
  fMCParticleName         (Names::gkMCPartBrn),
  fMCEventInfoName        (Names::gkMCEvtInfoBrn),
  fPileUpName             (Names::gkPileupInfoBrn),  
  fPFSuperClusterName     ("PFSuperClusters"),
  fPFMetName              ("PFMet"),
  fPFJetName              (Names::gkPFJetBrn),
  funcorrPFJetName        ("AKt5PFJets"),
  fGenJetName             ("AKT5GenJets"),
  fLeptonTagElectronsName ("HggLeptonTagElectrons"),
  fLeptonTagMuonsName     ("HggLeptonTagMuons"),
  fLeptonTagSoftElectronsName ("HggLeptonTagSoftElectrons"),
  fLeptonTagSoftMuonsName     ("HggLeptonTagSoftMuons"),

  fIsData                 (false),
  fIsCutBased             (false),
  fPhotonsFromBranch      (kTRUE),  
  fPVFromBranch           (kTRUE),
  fGoodElectronsFromBranch(kTRUE),
  fPFJetsFromBranch       (kTRUE),
  // ----------------------------------------
  // flag for synchronization, adds vertex variables
  // should be on for synching trees
  fDoSynching             (kFALSE),

  // ----------------------------------------
  // collections....
  fPhotons                (0),
  fPFPhotons              (0),
  fElectrons              (0),
  fConversions            (0),
  fPFConversions          (0),
  fTracks                 (0),
  fPileUpDen              (0),
  fPV                     (0),
  fBeamspot               (0),
  fPFCands                (0),
  fMCParticles            (0),
  fMCEventInfo            (0),
  fPileUp                 (0),
  fPFSuperClusters        (0),
  fPFJets                 (0),
  fGenJets                (0),
  funcorrPFJets           (0),

  fLeptonTagElectrons     (0),
  fLeptonTagMuons         (0),
  fPFNoPileUpCands        (0),
  fPFPileUpCands          (0),

  fLoopOnGoodElectrons    (kFALSE),
  fApplyElectronVeto      (kTRUE),  
  fWriteDiphotonTree      (kTRUE),
  fWriteSingleTree        (kTRUE),
  fEnablePFPhotons        (kTRUE),
  fExcludeSinglePrompt    (kFALSE),
  fExcludeDoublePrompt    (kFALSE),
  fEnableJets             (kFALSE),
  fEnableGenJets          (kFALSE),
  fApplyJetId             (kFALSE),
  fApplyLeptonTag         (kFALSE),
  fApplyVHLepTag          (kFALSE),
  fApplyVHHadTag          (kFALSE),
  fApplyVBFTag            (kFALSE),
  fApplyTTHTag            (kFALSE),
  fApplyBTag              (kFALSE),
  fApplyPFMetCorrections  (kFALSE),
  fFillClusterArrays      (kFALSE),
  fFillVertexTree         (kFALSE),
  fDo2012LepTag           (kFALSE),
  fVerbosityLevel         (0),
  fPhFixDataFile          (gSystem->Getenv("CMSSW_BASE") +
		           TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat")),
  fBeamspotWidth          (5.8),
  fTmpFile                (0),
  
  // JV: moved up the initializtion of fTupleName to avoid compilation warning
  fTupleName              ("hPhotonTree"),

  fElectronIDMVA(0),
  fElectronMVAWeights_Subdet0Pt10To20(""),
  fElectronMVAWeights_Subdet1Pt10To20(""),
  fElectronMVAWeights_Subdet2Pt10To20(""),
  fElectronMVAWeights_Subdet0Pt20ToInf(""),
  fElectronMVAWeights_Subdet1Pt20ToInf(""),
  fElectronMVAWeights_Subdet2Pt20ToInf(""),

  fTheRhoType(RhoUtilities::DEFAULT),

  fdor9rescale                   (false),
  fp0b                           (0.),
  fp1b                           (1.),
  fp0e                           (0.),
  fp1e                           (1.),

  fProcessedEvents(0)

{
  // Constructor
}

PhotonTreeWriter::~PhotonTreeWriter()
{
  // Destructor
  // Deal with the temporary file here?
  // fTmpFile->Write();
   fTmpFile->Close();
//   TString shellcmd = TString("rm ") + TString(fTmpFile->GetName());
//   delete fTmpFile;
//   cout << shellcmd.Data() << endl;
//   gSystem->Exec(shellcmd.Data());
  
}

//--------------------------------------------------------------------------------------------------
void PhotonTreeWriter::SlaveTerminate()
{
 
  if (hCiCTuple) fTmpFile->WriteTObject(hCiCTuple,hCiCTuple->GetName());
  if (hCiCTupleSingle) fTmpFile->WriteTObject(hCiCTuple,hCiCTupleSingle->GetName());
  
}


//--------------------------------------------------------------------------------------------------
void PhotonTreeWriter::Process()
{

  fProcessedEvents++;

  if (fVerbosityLevel > 0) {
    PhotonTreeWriter::LogEventInfo();
  }

  //if(GetEventHeader()->EvtNum()==9008 || GetEventHeader()->EvtNum()==9008 || GetEventHeader()->EvtNum()==9010){
  //  printf("check photontreewriter 0\n");
  // }
  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fGoodElectronName,   fGoodElectrons);

  // lepton tag collections
  if( fApplyLeptonTag || fApplyVHLepTag ) {
    LoadEventObject(fLeptonTagElectronsName, fLeptonTagElectrons);
    LoadEventObject(fLeptonTagMuonsName,     fLeptonTagMuons);
  }
  
  if( fApplyVHLepTag ) {
    LoadEventObject(fLeptonTagSoftElectronsName, fLeptonTagSoftElectrons);
    LoadEventObject(fLeptonTagSoftMuonsName,     fLeptonTagSoftMuons);
  }
  
  const BaseCollection *egcol = 0;
  if (fLoopOnGoodElectrons)
    egcol = fGoodElectrons;
  else
    egcol = fPhotons;
  if (egcol->GetEntries()<1)
    return;
  
  if (fEnablePFPhotons) LoadEventObject(fPFPhotonName,   fPFPhotons);
  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fConversionName,     fConversions);
  if ( fDoSynching ) LoadEventObject(fPFConversionName,     fPFConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);    
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadEventObject(fPFSuperClusterName, fPFSuperClusters);
  LoadEventObject(fPFMetName,          fPFMet);  
//   LoadEventObject(fPFNoPileUpName,     fPFNoPileUpCands);
//   LoadEventObject(fPFPileUpName,     fPFPileUpCands);

  if (fEnableJets){
    LoadEventObject(fPFJetName,        fPFJets);  
    //LoadEventObject(funcorrPFJetName,  funcorrPFJets);
    LoadBranch(funcorrPFJetName);
    //   if(!fIsData) LoadEventObject(fGenJetName,        fGenJets);
  }

  
  // ------------------------------------------------------------  
  // load event based information
  Int_t _numPU      = -1.;        // some sensible default values....
  Int_t _numPUminus = -1.;        // some sensible default values....
  Int_t _numPUplus  = -1.;        // some sensible default values....

  Float_t rho  = -99.;
  if( fPileUpDen->GetEntries() > 0 )
    rho  = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();//m
  
  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));//ming:what's this?
    
  if( !fIsData ) {
    LoadBranch(fMCParticleName);
    LoadBranch(fMCEventInfoName);
    LoadBranch(fPileUpName);
    if (fEnableGenJets) LoadEventObject(fGenJetName,        fGenJets);
  }  else fGenJets = NULL;
  
  if( !fIsData ) {
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing()==0) _numPU = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() == -1) _numPUminus = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() ==  1) _numPUplus  = puinfo->GetPU_NumInteractions();
    }
  }

  // in case of a MC event, try to find Higgs and Higgs decay Z poisition
  Float_t _pth     = -100.;
  Float_t _decayZ  = -100.;
  Float_t _genmass = -100.;
  if (!fIsData)
    FindHiggsPtAndZ(_pth, _decayZ, _genmass);

  
  fDiphotonEvent->genz = -999.;
  if (!fIsData) {
    fDiphotonEvent->genz = fMCParticles->At(0)->DecayVertex().Z();
  }
  
  fDiphotonEvent->mcprocid = -999;
  if (!fIsData) {
    fDiphotonEvent->mcprocid = fMCEventInfo->ProcessId();
    fDiphotonEvent->mcweight = fMCEventInfo->Weight();
  }

  Double_t _evt = GetEventHeader()->EvtNum();

  Double_t _spfMet = fPFMet->At(0)->SumEt();

  fDiphotonEvent->leptonTag = -1; // disabled
  fDiphotonEvent->VHLepTag = -1; // disabled
  fDiphotonEvent->VHHadTag = -1; // disabled

  // ====================================================
  // Vtx synching stuff...
  fDiphotonEvent->vtxInd1 = -1;
  fDiphotonEvent->vtxInd2 = -1;
  fDiphotonEvent->vtxInd3 = -1;

  fDiphotonEvent->vtxBestPtbal  = -1.;
  fDiphotonEvent->vtxBestPtasym = -1.;
  fDiphotonEvent->vtxBestSumpt2 = -1.;
  fDiphotonEvent->vtxBestP2Conv = -1.;

  fDiphotonEvent->vtxMva1 = -1.;
  fDiphotonEvent->vtxMva2 = -1.;
  fDiphotonEvent->vtxMva3 = -1.;

  fDiphotonEvent->vtxNleg1 = -1;
  fDiphotonEvent->vtxNleg2 = -1;
  fDiphotonEvent->vtxNconv = -1;
  // ====================================================


  fDiphotonEvent->rho = fPileUpDen->At(0)->RhoKt6PFJets();
  fDiphotonEvent->rho25 = fPileUpDen->At(0)->RhoRandomLowEta();
  fDiphotonEvent->rhoold = fPileUpDen->At(0)->Rho();
  fDiphotonEvent->genHiggspt = _pth;
  fDiphotonEvent->genHiggsZ = _decayZ;
  fDiphotonEvent->genmass = _genmass;  
  fDiphotonEvent->gencostheta = -99.;
  fDiphotonEvent->nVtx = fPV->GetEntries();
  fDiphotonEvent->bsX = fBeamspot->At(0)->X();
  fDiphotonEvent->bsY = fBeamspot->At(0)->Y();
  fDiphotonEvent->bsZ = fBeamspot->At(0)->Z();
  fDiphotonEvent->bsSigmaZ = fBeamspot->At(0)->SigmaZ();
  fDiphotonEvent->vtxX = (fDiphotonEvent->nVtx>0) ? fPV->At(0)->X() : -99.;
  fDiphotonEvent->vtxY = (fDiphotonEvent->nVtx>0) ? fPV->At(0)->Y() : -99.;  
  fDiphotonEvent->vtxZ = (fDiphotonEvent->nVtx>0) ? fPV->At(0)->Z() : -99.;
  fDiphotonEvent->numPU = _numPU;
  fDiphotonEvent->numPUminus = _numPUminus;
  fDiphotonEvent->numPUplus = _numPUplus;
  fDiphotonEvent->mass = -99.;
  fDiphotonEvent->ptgg = -99.;
  fDiphotonEvent->costheta = -99.;
  fDiphotonEvent->mt = -99.;
  fDiphotonEvent->cosphimet = -99.;
  fDiphotonEvent->mtele = -99.;
  fDiphotonEvent->cosphimetele = -99.;
  fDiphotonEvent->evt = GetEventHeader()->EvtNum();
  fDiphotonEvent->run = GetEventHeader()->RunNum();
  fDiphotonEvent->lumi = GetEventHeader()->LumiSec();
  fDiphotonEvent->evtcat = -99;
  fDiphotonEvent->nobj = fPhotons->GetEntries();
  fDiphotonEvent->pfmet = fPFMet->At(0)->Pt();
  fDiphotonEvent->pfmetphi = fPFMet->At(0)->Phi();
  fDiphotonEvent->pfmetx = fPFMet->At(0)->Px();
  fDiphotonEvent->pfmety = fPFMet->At(0)->Py();
  fDiphotonEvent->spfMet = _spfMet;
  fDiphotonEvent->masscor = -99.;
  fDiphotonEvent->masscorerr = -99.;
  fDiphotonEvent->masscorele = -99.;
  fDiphotonEvent->masscoreleerr = -99.;
  fDiphotonEvent->ismc = GetEventHeader()->IsMC();
  
  //jets
  const Jet *jet1 = 0;
  const Jet *jet2 = 0;
  const Jet *jetcentral = 0;

  fDiphotonEvent->jetleadNoIDpt   = -99.;
  fDiphotonEvent->jetleadNoIDeta  = -99.;
  fDiphotonEvent->jetleadNoIDphi  = -99.;
  fDiphotonEvent->jetleadNoIDmass = -99.;
  fDiphotonEvent->jet1pt   = -99.;
  fDiphotonEvent->jet1eta  = -99.;
  fDiphotonEvent->jet1phi  = -99.;
  fDiphotonEvent->jet1mass = -99.;
  fDiphotonEvent->jet2pt   = -99.;
  fDiphotonEvent->jet2eta  = -99.;
  fDiphotonEvent->jet2phi  = -99.;
  fDiphotonEvent->jet2mass = -99.;
  fDiphotonEvent->jetcentralpt   = -99.;
  fDiphotonEvent->jetcentraleta  = -99.;
  fDiphotonEvent->jetcentralphi  = -99.;
  fDiphotonEvent->jetcentralmass = -99.;
  fDiphotonEvent->dijetpt = -99.;
  fDiphotonEvent->dijeteta = -99.;
  fDiphotonEvent->dijetphi = -99.;
  fDiphotonEvent->dijetmass = -99.;   
  fDiphotonEvent->jetetaplus = -99.;
  fDiphotonEvent->jetetaminus = -99.;
  fDiphotonEvent->zeppenfeld = -99.;
  fDiphotonEvent->dphidijetgg = -99.;
  
  //uncorrected jets
//   const Jet *uncorrjet1 = 0;
//   const Jet *uncorrjet2 = 0;
//   const Jet *uncorrjetcentral = 0;

//   fDiphotonEvent->uncorrjet1pt   = -99.;
//   fDiphotonEvent->uncorrjet1eta  = -99.;
//   fDiphotonEvent->uncorrjet1phi  = -99.;
//   fDiphotonEvent->uncorrjet1mass = -99.;
//   fDiphotonEvent->uncorrjet2pt   = -99.;
//   fDiphotonEvent->uncorrjet2eta  = -99.;
//   fDiphotonEvent->uncorrjet2phi  = -99.;
//   fDiphotonEvent->uncorrjet2mass = -99.;
//   fDiphotonEvent->uncorrjetcentralpt   = -99.;
//   fDiphotonEvent->uncorrjetcentraleta  = -99.;
//   fDiphotonEvent->uncorrjetcentralphi  = -99.;
//   fDiphotonEvent->uncorrjetcentralmass = -99.;
//   fDiphotonEvent->diuncorrjetpt = -99.;
//   fDiphotonEvent->diuncorrjeteta = -99.;
//   fDiphotonEvent->diuncorrjetphi = -99.;
//   fDiphotonEvent->diuncorrjetmass = -99.; 


 
  Int_t nhitsbeforevtxmax = 1;
  if (!fApplyElectronVeto)
    nhitsbeforevtxmax = 999;  
  
  if (egcol->GetEntries()>=2) {
    const Particle *p1 = 0;
    const Particle *p2 = 0;
    const Photon *phHard = 0;
    const Photon *phSoft = 0;
    const Electron *ele1 = 0;
    const Electron *ele2 = 0;
    const SuperCluster *sc1 = 0;
    const SuperCluster *sc2 = 0;
   
    if (fLoopOnGoodElectrons) {
      ele1 = fGoodElectrons->At(0);
      ele2 = fGoodElectrons->At(1);
      p1 = ele1;
      p2 = ele2;
      sc1 = ele1->SCluster();
      sc2 = ele2->SCluster();
      phHard = PhotonTools::MatchedPhoton(ele1,fPhotons);
      phSoft = PhotonTools::MatchedPhoton(ele2,fPhotons);
    }
    else {
      phHard = fPhotons->At(0);
      phSoft = fPhotons->At(1);
      p1 = phHard;
      p2 = phSoft;
      sc1 = phHard->SCluster();
      sc2 = phSoft->SCluster();      
      ele1 = PhotonTools::MatchedElectron(phHard,fGoodElectrons);
      ele2 = PhotonTools::MatchedElectron(phSoft,fGoodElectrons);
    }
    
    const DecayParticle *conv1 = PhotonTools::MatchedConversion(sc1,fConversions,bsp,
								nhitsbeforevtxmax);
    const DecayParticle *conv2 = PhotonTools::MatchedConversion(sc2,fConversions,bsp,
								nhitsbeforevtxmax);
    
    const SuperCluster *pfsc1 = PhotonTools::MatchedPFSC(sc1,fPFPhotons, fElectrons);
    const SuperCluster *pfsc2 = PhotonTools::MatchedPFSC(sc2,fPFPhotons, fElectrons);
    
    const MCParticle *phgen1 = 0;
    const MCParticle *phgen2 = 0;
    if (!fIsData) {
      phgen1 = PhotonTools::MatchMC(p1,fMCParticles,!fApplyElectronVeto);
      phgen2 = PhotonTools::MatchMC(p2,fMCParticles,!fApplyElectronVeto);
    }
  
    if (fExcludeSinglePrompt && (phgen1 || phgen2))
      return;
    if (fExcludeDoublePrompt && (phgen1 && phgen2))
      return;
    
    if (!fLoopOnGoodElectrons && phHard->HasPV()) {
      fDiphotonEvent->vtxX = phHard->PV()->X();
      fDiphotonEvent->vtxY = phHard->PV()->Y();
      fDiphotonEvent->vtxZ = phHard->PV()->Z();
      fDiphotonEvent->vtxprob = phHard->VtxProb();
    }
    
    // fill Btag information... set to true if One jet fullfills 
    fDiphotonEvent -> btagJet1       = -1.;
    fDiphotonEvent -> btagJet1Pt     = -99.;
    fDiphotonEvent -> btagJet1Eta    = -99.;

    fDiphotonEvent -> btagJet2       = -1.;
    fDiphotonEvent -> btagJet2Pt     = -99.;
    fDiphotonEvent -> btagJet2Eta    = -99.;
   
    if( fApplyBTag && fEnableJets ) {
      float highTag     = 0.;
      float highJetPt   = 0.;
      float highJetEta  = 0.;

      float trailTag     = 0.;
      float trailJetPt   = 0.;
      float trailJetEta  = 0.;

      for (UInt_t ijet=0; ijet<fPFJets->GetEntries();++ijet) {
	const Jet *jet = fPFJets->At(ijet);
	if ( jet->Pt() < 20. || jet->AbsEta() > 2.4 ) continue;
	if ( jet->CombinedSecondaryVertexBJetTagsDisc() > highTag ) {
	  
	  trailTag    = highTag;
	  trailJetPt  = highJetPt;
	  trailJetEta = highJetEta;

	  highTag    = jet->CombinedSecondaryVertexBJetTagsDisc();
	  highJetPt  = jet->Pt();
	  highJetEta = jet->Eta();

	} else if ( jet->CombinedSecondaryVertexBJetTagsDisc() > trailTag ) {

	  trailTag    = jet->CombinedSecondaryVertexBJetTagsDisc();
	  trailJetPt  = jet->Pt();
	  trailJetEta = jet->Eta();
	}
      }
      fDiphotonEvent -> btagJet1    = highTag;
      fDiphotonEvent -> btagJet1Pt  = highJetPt;
      fDiphotonEvent -> btagJet1Eta = highJetEta;      

      fDiphotonEvent -> btagJet2    = trailTag;
      fDiphotonEvent -> btagJet2Pt  = trailJetPt;
      fDiphotonEvent -> btagJet2Eta = trailJetEta;      

    }

       
    //fill jet variables
    const Vertex *selvtx = fPV->At(0);
    if (!fLoopOnGoodElectrons && phHard->HasPV()) selvtx = phHard->PV();
    if (fEnableJets) {
      for (UInt_t ijet=0; ijet<fPFJets->GetEntries();++ijet) {
        const Jet *jet = fPFJets->At(ijet);
        if (jet->AbsEta()<4.7 && MathUtils::DeltaR(jet,p1)>0.5 && MathUtils::DeltaR(jet,p2)>0.5) {
          const PFJet *pfjet = dynamic_cast<const PFJet*>(jet);
          if (!pfjet) continue;
          //if (!JetTools::passPFLooseId(pfjet)) continue;
          if (fApplyJetId && !fJetId.passCut(pfjet,selvtx,fPV)) continue;
          if (!jet1) jet1 = jet;
          else if (!jet2) jet2 = jet;
          else if (!jetcentral && 0) jetcentral = jet;
        }
        if (jet1&&jet2&&jetcentral) break;
      }
      
      for(UInt_t ijet=0; ijet<fPFJets->GetEntries();++ijet){

	const Jet *jet = fPFJets->At(ijet);
	if (MathUtils::DeltaR(jet,p1)>0.5 && MathUtils::DeltaR(jet,p2)>0.5) {
	  fDiphotonEvent->jetleadNoIDpt   = jet->Pt();
	  fDiphotonEvent->jetleadNoIDeta  = jet->Eta();
	  fDiphotonEvent->jetleadNoIDphi  = jet->Phi();
	  fDiphotonEvent->jetleadNoIDmass = jet->Mass();
	  break;
	}
      }
      
      if (jet1) {
	fDiphotonEvent->jet1pt   = jet1->Pt();
	fDiphotonEvent->jet1eta  = jet1->Eta();
	fDiphotonEvent->jet1phi  = jet1->Phi();
	fDiphotonEvent->jet1mass = jet1->Mass();
      }
      
      if (jet2) {
	fDiphotonEvent->jet2pt   = jet2->Pt();
	fDiphotonEvent->jet2eta  = jet2->Eta();
	fDiphotonEvent->jet2phi  = jet2->Phi();
	fDiphotonEvent->jet2mass = jet2->Mass();
      }
      
      if (jetcentral) {
	fDiphotonEvent->jetcentralpt   = jetcentral->Pt();
	fDiphotonEvent->jetcentraleta  = jetcentral->Eta();
	fDiphotonEvent->jetcentralphi  = jetcentral->Phi();
	fDiphotonEvent->jetcentralmass = jetcentral->Mass();
      }
      
      if (jet1&&jet2){
	FourVectorM momjj = (jet1->Mom() + jet2->Mom());
	
	fDiphotonEvent->dijetpt =  momjj.Pt();
	fDiphotonEvent->dijeteta = momjj.Eta();
	fDiphotonEvent->dijetphi = momjj.Phi();
	fDiphotonEvent->dijetmass = momjj.M();    
	
	if (jet1->Eta()>jet2->Eta()) {
	  fDiphotonEvent->jetetaplus = jet1->Eta();
	  fDiphotonEvent->jetetaminus = jet2->Eta();
	}
	else {
	  fDiphotonEvent->jetetaplus = jet2->Eta();
	  fDiphotonEvent->jetetaminus = jet1->Eta();      
	}
      }
    }
   
    //added gen. info of whether a lep. or nutrino is from W or Z --Heng 02/14/2012 12:30 EST
    Double_t _fromZ = -99;
    Double_t _fromW = -99;
    Float_t _zpt = -99;
    Float_t _zEta = -99;
    Float_t _allZpt = -99;
    Float_t _allZEta = -99;
    
    if( !fIsData ){
      
      // loop over all GEN particles and look for nutrinos whoes mother is Z
      for(UInt_t j=0; j<fMCParticles->GetEntries(); ++j) {
	const MCParticle* p = fMCParticles->At(j);
	if( p->AbsPdgId()==23 ||p->AbsPdgId()==32 || p->AbsPdgId()==33 ) {
	  _allZpt=p->Pt();
	  _allZEta=p->Eta();
	  if (p->HasDaughter(12,kFALSE) || p->HasDaughter(14,kFALSE) || p->HasDaughter(16,kFALSE) ||p->HasDaughter(18,kFALSE) ) {
	  _fromZ=1;
	  _zpt=p->Pt();
	  _zEta=p->Eta();
	  }
	}
	else _fromW=1;
      }
    }
   
    fDiphotonEvent->fromZ = _fromZ;
    fDiphotonEvent->fromW = _fromW;
    fDiphotonEvent->zpt = _zpt;
    fDiphotonEvent->zEta = _zEta;
    fDiphotonEvent->allZpt = _allZpt;
    fDiphotonEvent->allZEta = _allZEta;

    Double_t _dphiMetgg = -99;
    Double_t _cosdphiMetgg = -99;
    Double_t _dphiPhPh = -99;

    Double_t _mass = -99.;
    Double_t _masserr = -99.;
    Double_t _masserrsmeared = -99.;
    Double_t _masserrwrongvtx = -99.;
    Double_t _masserrsmearedwrongvtx = -99.;
    Double_t _ptgg = -99.;
    Double_t _etagg = -99.;
    Double_t _phigg = -99.;    
    Double_t _costheta = -99.;
    PhotonTools::DiphotonR9EtaPtCats _evtcat = PhotonTools::kOctCat0;

    const Vertex* realVtx = NULL;

    if (phHard && phSoft) {
      _dphiMetgg = MathUtils::DeltaPhi((phHard->Mom()+phSoft->Mom()).Phi(),fPFMet->At(0)->Phi());
      _cosdphiMetgg = TMath::Cos(_dphiMetgg);
      _dphiPhPh = MathUtils::DeltaPhi((phHard->Mom()).Phi(),(phSoft->Mom()).Phi());
      _mass = (phHard->Mom()+phSoft->Mom()).M();
      _masserr = 0.5*TMath::Sqrt(phHard->EnergyErr()*phHard->EnergyErr() + phSoft->EnergyErr()*phSoft->EnergyErr());
      _masserrsmeared = 0.5*TMath::Sqrt(phHard->EnergyErrSmeared()*phHard->EnergyErrSmeared() + phSoft->EnergyErrSmeared()*phSoft->EnergyErrSmeared());
      _ptgg = (phHard->Mom()+phSoft->Mom()).Pt();
      _etagg = (phHard->Mom()+phSoft->Mom()).Eta();
      _phigg = (phHard->Mom()+phSoft->Mom()).Phi();
      _costheta = ThreeVector(phHard->Mom()).Unit().Dot(ThreeVector(phSoft->Mom()).Unit());
      _evtcat = PhotonTools::DiphotonR9EtaPtCat(phHard,phSoft);
      
      const Double_t dz = sqrt(2.0)*fBeamspotWidth;
      Double_t deltamvtx = VertexTools::DeltaMassVtx(phHard->CaloPos().X(),
							   phHard->CaloPos().Y(),
							   phHard->CaloPos().Z(),
							   phSoft->CaloPos().X(),
							   phSoft->CaloPos().Y(),
							   phSoft->CaloPos().Z(),
							   fDiphotonEvent->vtxX,
							   fDiphotonEvent->vtxY,
							   fDiphotonEvent->vtxZ,
							   dz);
      fDiphotonEvent->deltamvtx = deltamvtx;
      
      _masserrwrongvtx        = TMath::Sqrt(_masserr*_masserr + deltamvtx*deltamvtx);
      _masserrsmearedwrongvtx = TMath::Sqrt(_masserrsmeared*_masserrsmeared + deltamvtx*deltamvtx);
      


      // =================================================================================
      // this is for synching the Vtx stuff
      if (  fDoSynching ) {
	//fill conversion collection for vertex selection, adding single leg conversions if needed
	//note that momentum of single leg conversions needs to be recomputed from the track
	//as it is not filled properly
	DecayParticleOArr vtxconversions;
	if ( true ) {
	  vtxconversions.SetOwner(kTRUE);
	  for (UInt_t iconv=0; iconv<fConversions->GetEntries(); ++iconv) {
	    DecayParticle *conv = new DecayParticle(*fConversions->At(iconv));
	    vtxconversions.AddOwned(conv);
	  }
	  
	  for (UInt_t iconv=0; iconv<fPFConversions->GetEntries(); ++iconv) {
	    const DecayParticle *c = fPFConversions->At(iconv);
	    if (c->NDaughters()!=1) continue;
	    
	    DecayParticle *conv = new DecayParticle(*c);
	    const Track *trk = static_cast<const StableParticle*>(conv->Daughter(0))->Trk();
	    conv->SetMom(trk->Px(), trk->Py(), trk->Pz(), trk->P());
	    vtxconversions.AddOwned(conv);
	  }    
	}
	else {
	  for (UInt_t iconv=0; iconv<fConversions->GetEntries(); ++iconv) {
	    const DecayParticle *c = fConversions->At(iconv);
	    vtxconversions.Add(c);
	  }
	}
	
	
	const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));
	double tmpVtxProb = 0.;
	
	std::vector<int>    debugInds;  // this hold the Vtx indices for the first three guys
	std::vector<double> debugVals;  // this holds hte mva input/output for the best ranked
	std::vector<int>    debugConv;  // this holds number of legs for frist and second photon (0 if no conversion) and total number of conversions
	
		
	realVtx = fVtxTools.findVtxBasicRanking(phHard, phSoft, bsp, fPV,
						&vtxconversions,kTRUE,tmpVtxProb,
						&debugInds, &debugVals, &debugConv
						);
	
	fDiphotonEvent->vtxInd1 = debugInds[0];
	fDiphotonEvent->vtxInd2 = debugInds[1];
	fDiphotonEvent->vtxInd3 = debugInds[2];
	
	fDiphotonEvent->vtxConv1Z    = debugVals[0];
	fDiphotonEvent->vtxConv1DZ   = debugVals[1];
	fDiphotonEvent->vtxConv1Prob = debugVals[2];

	fDiphotonEvent->vtxConv2Z    = debugVals[3];
	fDiphotonEvent->vtxConv2DZ   = debugVals[4];
	fDiphotonEvent->vtxConv2Prob = debugVals[5];


	fDiphotonEvent->vtxBestPtbal  = debugVals[6];
	fDiphotonEvent->vtxBestPtasym = debugVals[7];
	fDiphotonEvent->vtxBestSumpt2 = debugVals[8];
	fDiphotonEvent->vtxBestP2Conv = debugVals[9];
	
	fDiphotonEvent->vtxMva1Z      = debugVals[10];
	fDiphotonEvent->vtxMva2Z      = debugVals[11];
	fDiphotonEvent->vtxMva3Z      = debugVals[12];

	fDiphotonEvent->vtxMva1 = debugVals[13];
	fDiphotonEvent->vtxMva2 = debugVals[14];
	fDiphotonEvent->vtxMva3 = debugVals[15];
	
	fDiphotonEvent->vtxNleg1    = debugConv[0];
	fDiphotonEvent->vtxNleg2    = debugConv[1];
	fDiphotonEvent->vtxConvIdx1 = debugConv[2];
	fDiphotonEvent->vtxConvIdx2 = debugConv[3];
	
	fDiphotonEvent->vtxNconv = debugConv[4];

	if( false ) {
	  printf("---------------------------------------------------------------------\n");
	  printf("FINAL\nvtx %i: ptbal = %5f, ptasym = %5f, logsumpt2 = %5f, limpulltoconv = %5f, nconv = %i, mva = %5f\n",fDiphotonEvent->vtxInd1,fDiphotonEvent->vtxBestPtbal,fDiphotonEvent->vtxBestPtasym,fDiphotonEvent->vtxBestSumpt2,fDiphotonEvent->vtxBestP2Conv,fDiphotonEvent->vtxNconv,fDiphotonEvent->vtxMva1);
	  printf("---------------------------------------------------------------------\n");
	}
      }
      
      // end of Vtx synching stuff ...
      // =================================================================================
      	
	PFJetOArr pfjets;
	if (fEnableJets){
	  if (jet1 && jet2) {
	    fDiphotonEvent->zeppenfeld = TMath::Abs(_etagg - 0.5*(jet1->Eta()+jet2->Eta()));
	    fDiphotonEvent->dphidijetgg = MathUtils::DeltaPhi( (jet1->Mom()+jet2->Mom()).Phi(), _phigg );
	  }
	  
	  //PFJetOArr pfjets;
	  for (UInt_t ijet=0; ijet<fPFJets->GetEntries(); ++ijet) {
	    const PFJet *pfjet = dynamic_cast<const PFJet*>(fPFJets->At(ijet));
	    if (pfjet && MathUtils::DeltaR(*pfjet,*phHard)>0.3 && MathUtils::DeltaR(*pfjet,*phSoft)>0.3) pfjets.Add(pfjet);
	  }
	}
	
	PFCandidateOArr pfcands;
	for (UInt_t icand=0; icand<fPFCands->GetEntries(); ++icand) {
	  const PFCandidate *pfcand = fPFCands->At(icand);
	  if (MathUtils::DeltaR(*pfcand,*phHard)>0.1 && MathUtils::DeltaR(*pfcand,*phSoft)>0.1) pfcands.Add(pfcand);
	}      
      
	const Vertex *firstvtx = fPV->At(0);
	const Vertex *selvtx = fPV->At(0);
	
	if (!fLoopOnGoodElectrons && phHard->HasPV()) {
	  selvtx = phHard->PV();
	}
	
	if (0) //disable for now for performance reasons
	  {
	    Met mmet = fMVAMet.GetMet(  false,
					0.,0.,0.,
					0.,0.,0.,
					fPFMet->At(0),
					&pfcands,selvtx,fPV, fPileUpDen->At(0)->Rho(),
					&pfjets,
					int(fPV->GetEntries()),
					kFALSE);      
	    
	    TMatrixD *metcov = fMVAMet.GetMetCovariance();
	    
	    ThreeVector fullmet(mmet.Px() - phHard->Px() - phSoft->Px(),
				mmet.Py() - phHard->Py() - phSoft->Py(),
				0.);
	    
	    ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > mcov;
	    mcov(0,0) = (*metcov)(0,0);
	    mcov(0,1) = (*metcov)(0,1);
	    mcov(1,0) = (*metcov)(1,0);
	    mcov(1,1) = (*metcov)(1,1);
	    ROOT::Math::SVector<double,2> vmet;
	    vmet(0) = fullmet.X();
	    vmet(1) = fullmet.Y();
	    mcov.Invert();
	    Double_t metsig = sqrt(ROOT::Math::Similarity(mcov,vmet));
	    
	    fDiphotonEvent->mvametsel = fullmet.Rho();
	    fDiphotonEvent->mvametselphi = fullmet.Phi();
	    fDiphotonEvent->mvametselx = fullmet.X();
	    fDiphotonEvent->mvametsely = fullmet.Y();
	    fDiphotonEvent->mvametselsig = metsig;
	  }
	
	if (0) //disable for now for performance reasons
	  {
	    Met mmet = fMVAMet.GetMet(  false,
					0.,0.,0.,
					0.,0.,0.,
					fPFMet->At(0),
					&pfcands,firstvtx,fPV, fPileUpDen->At(0)->Rho(),
					&pfjets,
					int(fPV->GetEntries()),
					kFALSE);      
	    
	    TMatrixD *metcov = fMVAMet.GetMetCovariance();
	    
	    ThreeVector fullmet(mmet.Px() - phHard->Px() - phSoft->Px(),
				mmet.Py() - phHard->Py() - phSoft->Py(),
				0.);
	    
	    ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > mcov;
	    mcov(0,0) = (*metcov)(0,0);
	    mcov(0,1) = (*metcov)(0,1);
	    mcov(1,0) = (*metcov)(1,0);
	    mcov(1,1) = (*metcov)(1,1);
	    ROOT::Math::SVector<double,2> vmet;
	    vmet(0) = fullmet.X();
	    vmet(1) = fullmet.Y();
	    mcov.Invert();
	    Double_t metsig = sqrt(ROOT::Math::Similarity(mcov,vmet));
	    
	    fDiphotonEvent->mvametfirst = fullmet.Rho();
	    fDiphotonEvent->mvametfirstphi = fullmet.Phi();
	    fDiphotonEvent->mvametfirstx = fullmet.X();
	    fDiphotonEvent->mvametfirsty = fullmet.Y();
	    fDiphotonEvent->mvametfirstsig = metsig;
	  }      
	
    }
    
    fDiphotonEvent->corrpfmet = -99.;
    fDiphotonEvent->corrpfmetphi = -99.;
    fDiphotonEvent->corrpfmetx = -99.;
    fDiphotonEvent->corrpfmety = -99.;
    
    Met *corrMet =NULL;
    
    if (fApplyPFMetCorrections){
      corrMet = new Met(fPFMet->At(0)->Px(),fPFMet->At(0)->Py());
      
      if (!fIsData){
	PFMetCorrectionTools::correctMet(corrMet,phHard,phSoft,1,0,funcorrPFJets,fGenJets,fPFJets,_evt);
	PFMetCorrectionTools::shiftMet(corrMet,fIsData,_spfMet);
      }
      else {
	PFMetCorrectionTools::shiftMet(corrMet,fIsData,_spfMet);
	PFMetCorrectionTools::correctMet(corrMet,phHard,phSoft,0,1,funcorrPFJets,fGenJets,fPFJets,_evt);
      }    
      
      fDiphotonEvent->corrpfmet = corrMet->Pt();
      fDiphotonEvent->corrpfmetphi = corrMet->Phi();
      fDiphotonEvent->corrpfmetx = corrMet->Px();
      fDiphotonEvent->corrpfmety = corrMet->Py();
      
      delete corrMet;
    }
    
    
    Float_t _massele = -99.;
    Float_t _ptee = -99.;
    Float_t _costhetaele = -99.;
    if (ele1 && ele2) {
      _massele = (ele1->Mom()+ele2->Mom()).M();
      _ptee = (ele1->Mom()+ele2->Mom()).Pt();
      _costhetaele = ThreeVector(ele1->Mom()).Unit().Dot(ThreeVector(ele2->Mom()).Unit());      
    }
    
    Float_t _gencostheta = -99.;
    if (phgen1 && phgen2) {
      _gencostheta = ThreeVector(phgen1->Mom()).Unit().Dot(ThreeVector(phgen2->Mom()).Unit());
    }
    
    fDiphotonEvent->gencostheta = _gencostheta;
    fDiphotonEvent->dphiMetgg = _dphiMetgg;
    fDiphotonEvent->cosdphiMetgg = _cosdphiMetgg;
    fDiphotonEvent->dphiPhPh = _dphiPhPh;
    fDiphotonEvent->mass = _mass;
    fDiphotonEvent->masserr = _masserr;
    fDiphotonEvent->masserrsmeared = _masserrsmeared;
    fDiphotonEvent->masserrwrongvtx = _masserrwrongvtx;
    fDiphotonEvent->masserrsmearedwrongvtx = _masserrsmearedwrongvtx;    
    fDiphotonEvent->ptgg = _ptgg;
    fDiphotonEvent->etagg = _etagg;
    fDiphotonEvent->phigg = _phigg;
    fDiphotonEvent->costheta =  _costheta;;
    fDiphotonEvent->massele = _massele;
    fDiphotonEvent->ptee = _ptee;
    fDiphotonEvent->costhetaele =  _costhetaele;    
    fDiphotonEvent->evtcat = _evtcat;
    
    fDiphotonEvent->photons[0].SetVars(phHard, conv1, ele1, pfsc1, phgen1,
                                       fPhfixph, fPhfixele, fTracks, fPV,
                                       fPFCands, rho, fFillClusterArrays,
                                       fPhotons, fPFSuperClusters,
                                       fElectrons, fConversions, bsp,fIsData, fdor9rescale, fp0b, fp1b,fp0e,fp1e,
                                       fApplyElectronVeto, realVtx);
    fDiphotonEvent->photons[1].SetVars(phSoft, conv2, ele2, pfsc2, phgen2,
                                       fPhfixph, fPhfixele, fTracks, fPV,
                                       fPFCands, rho, fFillClusterArrays,
                                       fPhotons, fPFSuperClusters,
                                       fElectrons, fConversions,
                                       bsp,fIsData, fdor9rescale, fp0b, fp1b,fp0e,fp1e, fApplyElectronVeto, realVtx);
    
    Float_t ph1ecor    = fDiphotonEvent->photons[0].Ecor();
    Float_t ph1ecorerr = fDiphotonEvent->photons[0].Ecorerr();
    Float_t ph2ecor    = fDiphotonEvent->photons[1].Ecor();
    Float_t ph2ecorerr = fDiphotonEvent->photons[1].Ecorerr();
    
    Float_t ph1ecorele    = fDiphotonEvent->photons[0].Ecorele();
    Float_t ph1ecoreleerr = fDiphotonEvent->photons[0].Ecoreleerr();
    Float_t ph2ecorele    = fDiphotonEvent->photons[1].Ecorele();
    Float_t ph2ecoreleerr = fDiphotonEvent->photons[1].Ecoreleerr();
    
    fDiphotonEvent->masscor    = TMath::Sqrt(2.0*ph1ecor*ph2ecor*(1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscorerr = 0.5*fDiphotonEvent->masscor*
      TMath::Sqrt(ph1ecorerr*ph1ecorerr/ph1ecor/ph1ecor + ph2ecorerr*ph2ecorerr/ph2ecor/ph2ecor);
    
    fDiphotonEvent->masscorele = TMath::Sqrt(2.0*ph1ecorele*ph2ecorele*
					     (1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscoreleerr = 0.5*fDiphotonEvent->masscorele*
      TMath::Sqrt(ph1ecoreleerr*ph1ecoreleerr/ph1ecorele/ph1ecorele +
		  ph2ecoreleerr*ph2ecoreleerr/ph2ecorele/ph2ecorele);    
    
    //printf("r9 = %5f, photon sigieie = %5f, seed sigieie = %5f\n",phHard->R9(),
    //       phHard->CoviEtaiEta(),sqrt(phHard->SCluster()->Seed()->CoviEtaiEta()));
    
    // MuonStuff
    fDiphotonEvent-> muonPt  = -99.;
    fDiphotonEvent-> muonEta  = -99.;
    fDiphotonEvent-> muonPhi  = -99.;
    fDiphotonEvent-> muDR1  = -99.;
    fDiphotonEvent-> muDR2  = -99.;
    fDiphotonEvent-> muIso1  = -99.;
    fDiphotonEvent-> muIso2  = -99.;
    fDiphotonEvent-> muIso3  = -99.;
    fDiphotonEvent-> muIso4  = -99.;
    fDiphotonEvent-> muD0  = -99.;
    fDiphotonEvent-> muDZ  = -99.;
    fDiphotonEvent-> muChi2  = -99.;   
    fDiphotonEvent-> muNhits = -99;
    fDiphotonEvent-> muNpixhits = -99;
    fDiphotonEvent-> muNegs = -99;
    fDiphotonEvent-> muNMatch = -99;

    // Dimuon Stuff
    fDiphotonEvent-> mu1Pt   = -99.;
    fDiphotonEvent-> mu1Eta  = -99.;
    fDiphotonEvent-> mu1Phi  = -99.;
    fDiphotonEvent-> mu2Pt   = -99.;
    fDiphotonEvent-> mu2Eta  = -99.;
    fDiphotonEvent-> mu2Phi  = -99.;
    
    // Electron Stuff
    fDiphotonEvent-> elePt = -99.;
    fDiphotonEvent-> eleEta = -99.;
    fDiphotonEvent-> elePhi = -99.;
    fDiphotonEvent-> eleSCEta = -99.;     
    fDiphotonEvent-> eleIso1 = -99.;
    fDiphotonEvent-> eleIso2 = -99.;
    fDiphotonEvent-> eleIso3 = -99.;
    fDiphotonEvent-> eleIso4 = -99.;
    fDiphotonEvent-> eleDist = -99.;
    fDiphotonEvent-> eleDcot = -99.;
    fDiphotonEvent-> eleCoviee = -99.;
    fDiphotonEvent-> eleDphiin = -99.;
    fDiphotonEvent-> eleDetain = -99.;
    fDiphotonEvent-> eleDR1 = -99.;
    fDiphotonEvent-> eleDR2 = -99.;
    fDiphotonEvent-> eleMass1 = -99.;
    fDiphotonEvent-> eleMass2 = -99.;
    fDiphotonEvent-> eleNinnerHits = -99;
    fDiphotonEvent-> eleIdMva = -99.;
    
    // Dielectron Stuff
    fDiphotonEvent-> ele1Pt   = -99.;
    fDiphotonEvent-> ele1Eta  = -99.;
    fDiphotonEvent-> ele1Phi  = -99.;
    fDiphotonEvent-> ele2Pt   = -99.;
    fDiphotonEvent-> ele2Eta  = -99.;
    fDiphotonEvent-> ele2Phi  = -99.;

    
    if( fApplyLeptonTag ) {
      ApplyLeptonTag(phHard, phSoft, selvtx);
    }

    if( fApplyVHLepTag ) {
      ApplyVHLepTag(phHard, phSoft, selvtx);
    }

    fDiphotonEvent->costhetastar = -99.;

    if( fApplyVHHadTag && jet1 && jet2) {
      ApplyVHHadTag(phHard, phSoft, selvtx, jet1, jet2);
    }

    //vbf tag
    fDiphotonEvent->vbfTag = -1;
    fDiphotonEvent->vbfbdt = -99;
    if(fApplyVBFTag){
      fDiphotonEvent->vbfTag = 0;
      jet1pt_vbf = fDiphotonEvent->jet1pt;
      jet2pt_vbf = fDiphotonEvent->jet2pt;
      deltajeteta_vbf = TMath::Abs(fDiphotonEvent->jet1eta - fDiphotonEvent->jet2eta);
      dijetmass_vbf = fDiphotonEvent->dijetmass;
      zeppenfeld_vbf = TMath::Abs(fDiphotonEvent->zeppenfeld);
      dphidijetgg_vbf = TMath::Abs(fDiphotonEvent->dphidijetgg);
      diphoptOverdiphomass_vbf = fDiphotonEvent->ptgg/fDiphotonEvent->mass;
      pho1ptOverdiphomass_vbf = phHard->Pt()/fDiphotonEvent->mass;
      pho2ptOverdiphomass_vbf = phSoft->Pt()/fDiphotonEvent->mass;
      
      if(jet1 && jet2 && (pho1ptOverdiphomass_vbf > 40./120.) && (pho2ptOverdiphomass_vbf > 30./120.) && (jet1pt_vbf > 30) && (jet2pt_vbf > 20) && (dijetmass_vbf > 250)){
	fDiphotonEvent->vbfbdt = fMVAVBF.GetMVAbdtValue(jet1pt_vbf,jet2pt_vbf,deltajeteta_vbf,dijetmass_vbf,zeppenfeld_vbf,dphidijetgg_vbf,diphoptOverdiphomass_vbf,pho1ptOverdiphomass_vbf,pho2ptOverdiphomass_vbf);
	
	if((fDiphotonEvent->vbfbdt>=0.985) && (fDiphotonEvent->vbfbdt<=1)){
	  fDiphotonEvent->vbfTag = 2;
	}
	if((fDiphotonEvent->vbfbdt>=0.93) && (fDiphotonEvent->vbfbdt<0.985)){
	  fDiphotonEvent->vbfTag = 1;
	}
      }
    } // End of VBF tag stuff
    
    // ttH tag stuff
    fDiphotonEvent->tthTag = -1;
    if (fApplyTTHTag && phHard && phSoft && selvtx) {
      ApplyTTHTag(phHard, phSoft, selvtx);
    }

    fDiphotonEvent->numJets  = -99;
    fDiphotonEvent->numBJets = -99;
    fDiphotonEvent->bjetcsv = -99;

    if (fDoSynching) {
      double minJetPt = 20.;
      double maxJetAbsEta = 4.7;
      fDiphotonEvent->numJets  = NumberOfJets (phHard, phSoft, selvtx, minJetPt,
                                               maxJetAbsEta);
      fDiphotonEvent->numBJets = NumberOfBJets(phHard, phSoft, selvtx, minJetPt,
                                               maxJetAbsEta);
      /// The lines below cause a segmentation violation in the call to 
      /// GetSelectedPFJets
//       DeltaRVetoVector vetos(2);
//       vetos.push_back(DeltaRVeto(static_cast<const Particle*>(phHard), 0.5));
//       vetos.push_back(DeltaRVeto(static_cast<const Particle*>(phSoft), 0.5));
//       // for (unsigned
//       PFJetVector *jets = GetSelectedPFJets(vetos, *selvtx, minJetPt, 
//                                             maxJetAbsEta);
//       fDiphotonEvent->numJets = jets->size();
//       delete jets;
      
    }

    
    //printf("vbfbdt:%f\n",fDiphotonEvent->vbfbdt);
    if (fWriteDiphotonTree)
      hCiCTuple->Fill();  
    
    if (fFillVertexTree && phHard && phSoft) {
      for (UInt_t ivtx = 0; ivtx<fPV->GetEntries(); ++ivtx) {
	const Vertex *v = fPV->At(ivtx);
	fDiphotonVtx->SetVars(v,phHard,phSoft,fPFCands,ivtx,fPV->GetEntries(),_decayZ);
	hVtxTree->Fill();
      }
    }
  }
 
  if (!fWriteSingleTree)
    return;
  
  for (UInt_t iph = 0; iph<egcol->GetEntries(); ++iph) {
    const Particle *p = 0;
    const Photon *ph = 0;
    const Electron *ele = 0;
    const SuperCluster *sc = 0;
    if (fLoopOnGoodElectrons) {
      ele = fGoodElectrons->At(iph);
      p = ele;
      sc = ele->SCluster();
      ph = PhotonTools::MatchedPhoton(ele,fPhotons);
    }
    else {
      ph = fPhotons->At(iph);
      p = ph;
      sc = ph->SCluster();
      ele = PhotonTools::MatchedElectron(ph,fGoodElectrons);    
    }
    
    const DecayParticle *conv = PhotonTools::MatchedConversion(sc,fConversions,bsp,
							       nhitsbeforevtxmax);
    const SuperCluster  *pfsc = PhotonTools::MatchedPFSC(sc,fPFPhotons, fElectrons);
        
    if (!fLoopOnGoodElectrons && ph->HasPV()) {
      fDiphotonEvent->vtxZ = ph->PV()->Z();
    }

    const MCParticle *phgen = 0;
    if( !fIsData ) {
      phgen = PhotonTools::MatchMC(p,fMCParticles,!fApplyElectronVeto);
    }

    if (fExcludeSinglePrompt && phgen) return;

    fDiphotonEvent->mt = -99.;
    fDiphotonEvent->cosphimet = -99.;
    fDiphotonEvent->mtele = -99.;
    fDiphotonEvent->cosphimetele = -99.;

    if (ph) {
      fDiphotonEvent->cosphimet = TMath::Cos(ph->Phi()-fPFMet->At(0)->Phi());
      fDiphotonEvent->mt = TMath::Sqrt(2.0*fPFMet->At(0)->Pt()*ph->Pt()*
				       (1.0-fDiphotonEvent->cosphimet));
    }
    
    if (ele) {
      fDiphotonEvent->cosphimetele = TMath::Cos(ele->Phi()-fPFMet->At(0)->Phi());
      fDiphotonEvent->mtele        = TMath::Sqrt(2.0*fPFMet->At(0)->Pt()*ele->Pt()*
						 (1.0-fDiphotonEvent->cosphimetele));      
    }
    
    fSinglePhoton->SetVars(ph, conv, ele, pfsc, phgen, fPhfixph, fPhfixele,
                           fTracks, fPV, fPFCands, rho, fFillClusterArrays,
                           fPhotons, fPFSuperClusters, fElectrons, fConversions,
                           bsp, fIsData, fdor9rescale, fp0b, fp1b,fp0e,fp1e, fApplyElectronVeto);
    hCiCTupleSingle->Fill();
  }



  return;
}

//--------------------------------------------------------------------------------------------------
void PhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  if( fApplyLeptonTag || fApplyVHLepTag ) {
    ReqEventObject(fLeptonTagElectronsName,    fLeptonTagElectrons,    false);  
    ReqEventObject(fLeptonTagMuonsName,        fLeptonTagMuons,        false);  
  }

  if( fApplyVHLepTag ) {
    ReqEventObject(fLeptonTagSoftElectronsName,    fLeptonTagSoftElectrons,    false);  
    ReqEventObject(fLeptonTagSoftMuonsName,        fLeptonTagSoftMuons,        false);  
  }

//   ReqEventObject(fPFNoPileUpName,     fPFNoPileUpCands,    false);
//   ReqEventObject(fPFPileUpName,     fPFPileUpCands,    false);

  ReqEventObject(fPhotonBranchName,fPhotons,      fPhotonsFromBranch);
  if (fEnablePFPhotons) ReqEventObject(fPFPhotonName,fPFPhotons,      true);
  ReqEventObject(fTrackBranchName, fTracks,       true);
  ReqEventObject(fElectronName,    fElectrons,    true);  
  ReqEventObject(fGoodElectronName,fGoodElectrons,fGoodElectronsFromBranch);  
  ReqEventObject(fPileUpDenName,   fPileUpDen,    true);
  ReqEventObject(fPVName,          fPV,           fPVFromBranch);
  ReqEventObject(fConversionName,  fConversions,  true);
  if ( fDoSynching ) ReqEventObject(fPFConversionName,     fPFConversions,  true);
  ReqEventObject(fBeamspotName,    fBeamspot,     true);
  ReqEventObject(fPFCandName,      fPFCands,      true);
  ReqEventObject(fPFSuperClusterName,fPFSuperClusters,true);
  ReqEventObject(fPFMetName,       fPFMet,        true);
  if (fEnableJets){
    ReqEventObject(fPFJetName,       fPFJets,       fPFJetsFromBranch);
    ReqBranch(funcorrPFJetName, funcorrPFJets);
    //   if (!fIsData) ReqEventObject(fGenJetName, fGenJets, true);
  }
  if (!fIsData) {
    ReqBranch(fPileUpName,         fPileUp);
    ReqBranch(fMCParticleName,     fMCParticles);
    ReqBranch(fMCEventInfoName,    fMCEventInfo);
    if (fEnableGenJets) ReqEventObject(fGenJetName, fGenJets, true);
  }
  if (fIsData) {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") +
      TString("/src/MitPhysics/data/PhotonFixGRPV22.dat");
  }
  else {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") +
      TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat");
  }

  fPhfixph.initialise("4_2",std::string(fPhFixDataFile));
  fPhfixele.initialise("4_2e",std::string(fPhFixDataFile));
  
//   fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
//                       TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2_42.root")))
//                       );  
  
  fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1cov_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2cov_52.root")))
                      );                      
                 
  fJetId.Initialize(JetIDMVA::kMedium,
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                          JetIDMVA::kCut,
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")))

                          );

  fMVAVBF.InitializeMVA();
                      

  if( fDoSynching ) {
    fVtxTools.InitP(2);
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               ElectronIDMVA::kIDEGamma2012NonTrigV1,
			       fTheRhoType);
  }


  fDiphotonEvent = new PhotonTreeWriterDiphotonEvent;
  fSinglePhoton  = new PhotonTreeWriterPhoton<16>;
  
  fTmpFile = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  
  if (fWriteDiphotonTree) {
    hCiCTuple = new TTree(fTupleName.Data(),fTupleName.Data());
    hCiCTuple->SetAutoSave(300e9);
    hCiCTuple->SetDirectory(fTmpFile);
  }
  TString singlename = fTupleName + TString("Single");
  if (fWriteSingleTree) {
    hCiCTupleSingle = new TTree(singlename,singlename);
    hCiCTupleSingle->SetAutoSave(300e9);
    hCiCTupleSingle->SetDirectory(fTmpFile);
  }
    
  //make flattish tree from classes so we don't have to rely on dictionaries for reading later
  TClass *eclass = TClass::GetClass("mithep::PhotonTreeWriterDiphotonEvent");
  TClass *pclass = TClass::GetClass("mithep::PhotonTreeWriterPhoton<16>");
  TList  *elist  = eclass->GetListOfDataMembers();
  TList  *plist  = pclass->GetListOfDataMembers();
    
  for (int i=0; i<elist->GetEntries(); ++i) {
    const TDataMember *tdm = static_cast<const TDataMember*>(elist->At(i));//ming
    if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
    TString typestring;
    if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
    else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
    else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
    else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
    else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
    else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
    else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
    else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
    else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
    else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
    else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
    else continue;
    //printf("%s %s: %i\n",tdm->GetTypeName(),tdm->GetName(),int(tdm->GetOffset()));
    Char_t *addr = (Char_t*)fDiphotonEvent;//ming:?
    assert(sizeof(Char_t)==1);
    if (fWriteDiphotonTree) hCiCTuple->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
    if (fWriteSingleTree) hCiCTupleSingle->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
  }

  for (int iph=0; iph<2; ++iph) {
    for (int i=0; i<plist->GetEntries(); ++i) {
      const TDataMember *tdm = static_cast<const TDataMember*>(plist->At(i));
      if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
      TString typestring;
      if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
      else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
      else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
      else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
      else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
      else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
      else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
      else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
      else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
      else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
      else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
      else continue;
      //printf("%s\n",tdm->GetTypeName());
      TString varname = TString::Format("ph%d.%s",iph+1,tdm->GetName());
      if (tdm->GetArrayDim()==1) {
        varname = TString::Format("%s[%i]",varname.Data(),tdm->GetMaxIndex(0));
      }
      
      //printf("typename = %s, arraydim = %i, arraysize = %i,varname = %s\n", tdm->GetTypeName(), tdm->GetArrayDim(), tdm->GetMaxIndex(0), varname.Data());
      
      Char_t *addr = (Char_t*)&fDiphotonEvent->photons[iph];
      assert(sizeof(Char_t)==1);
      if (fWriteDiphotonTree) hCiCTuple->Branch(varname,addr+tdm->GetOffset(),TString::Format("%s/%s",varname.Data(),typestring.Data()));
      
      if (iph==0) {
        TString singlename = TString::Format("ph.%s",tdm->GetName());
        if (tdm->GetArrayDim()==1) {
          singlename = TString::Format("%s[%i]",singlename.Data(),tdm->GetMaxIndex(0));
        }       
        Char_t *addrsingle = (Char_t*)fSinglePhoton;
        if (fWriteSingleTree) hCiCTupleSingle->Branch(singlename,addrsingle+tdm->GetOffset(),TString::Format("%s/%s",singlename.Data(),typestring.Data()));
      }
    }
  }
  
  if (fWriteDiphotonTree)
    AddOutput(hCiCTuple);
  if (fWriteSingleTree)
    AddOutput(hCiCTupleSingle);
  
  if (fFillVertexTree) {
    fDiphotonVtx = new PhotonTreeWriterVtx;
    hVtxTree = new TTree("hVtxTree","hVtxTree");
    hVtxTree->SetAutoSave(300e9);
    AddOutput(hVtxTree);
    
    TClass *vclass = TClass::GetClass("mithep::PhotonTreeWriterVtx");
    TList  *vlist  = vclass->GetListOfDataMembers();

    for (int i=0; i<vlist->GetEntries(); ++i) {
      const TDataMember *tdm = static_cast<const TDataMember*>(vlist->At(i));
      if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
      TString typestring;
      if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
      else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
      else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
      else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
      else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
      else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
      else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
      else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
      else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
      else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
      else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
      else continue;
      //printf("%s %s: %i\n",tdm->GetTypeName(),tdm->GetName(),int(tdm->GetOffset()));
      Char_t *addr = (Char_t*)fDiphotonVtx;
      assert(sizeof(Char_t)==1);
      hVtxTree->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
    }
  }
  
}

// ----------------------------------------------------------------------------------------
// some helpfer functions....
void PhotonTreeWriter::FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ, Float_t& mass)
{
  pt     = -999.;
  decayZ = -999.;
  mass   = -999.;

  // loop over all GEN particles and look for status 1 photons
  for(UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle* p = fMCParticles->At(i);
    if (p->Is(MCParticle::kH) || (!fApplyElectronVeto &&
				  (p->AbsPdgId()==23 || p->AbsPdgId()==24))) {
      pt     = p->Pt();
      decayZ = p->DecayVertex().Z();
      mass   = p->Mass();
      break;
    }
  }
  
  return;
 }


Float_t PhotonTreeWriter::GetEventCat(PhotonTools::CiCBaseLineCats cat1,
                                      PhotonTools::CiCBaseLineCats cat2) {
  
  bool ph1IsEB  = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat2);
  bool ph2IsEB  = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat2);

  bool ph1IsHR9 = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat3);
  bool ph2IsHR9 = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat3);
  
  if( ph1IsEB && ph2IsEB )
    return ( ph1IsHR9 && ph2IsHR9 ? 0. : 1.);
  
  return ( ph1IsHR9 && ph2IsHR9 ? 2. : 3.);
}

template <int NClus>
void
PhotonTreeWriterPhoton<NClus>::SetVars(const Photon *p,
                                       const DecayParticle *c,
                                       const Electron *ele,
                                       const SuperCluster *pfsc,
                                       const MCParticle *m,
                                       PhotonFix &phfixph,
                                       PhotonFix &phfixele,
                                       const TrackCol* trackCol,
                                       const VertexCol* vtxCol,
                                       const PFCandidateCol* fPFCands,
                                       Double_t rho,
                                       Bool_t fillclusterarrays,
                                       const PhotonCol* ps,
                                       const SuperClusterCol* scs,
                                       const ElectronCol* els,
                                       const DecayParticleCol *convs,
                                       const BaseVertex *bs,bool isdata,
				       bool dor9rescale, double p0b, double p1b, double  p0e,  double  p1e,
                                       Bool_t applyElectronVeto,
                                       const Vertex* realVtx)
{

  const SuperCluster *s = 0;
  if (p)
    s = p->SCluster();
  else
    s = ele->SCluster();
  const BasicCluster *b = s->Seed();
  const BasicCluster *b2 = 0;
  Double_t ebcmax = -99.;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() > ebcmax && bc !=b) {
      b2 = bc;
      ebcmax = bc->Energy();
    }
  }
  
  const BasicCluster *bclast = 0;
  Double_t ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b) {
      bclast = bc;
      ebcmin = bc->Energy();
    }
  }

  const BasicCluster *bclast2 = 0;
  ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b && bc!=bclast) {
      bclast2 = bc;
      ebcmin = bc->Energy();
    }
  }
  
  if (p) {
    hasphoton = kTRUE;
    index = PhotonTreeWriter::IndexOfElementInCollection(p, ps);
    
    e = p->E();
    pt = p->Pt();
    eta = p->Eta();
    phi = p->Phi();
    r9 = p->R9();
    e3x3 = p->E33();
    e5x5 = p->E55();
    hovere = p->HadOverEm();
    hoveretower = p->HadOverEmTow();
    sigietaieta = p->CoviEtaiEta();      
    phcat = PhotonTools::CiCBaseLineCat(p);
    eerr = p->EnergyErr();
    eerrsmeared = p->EnergyErrSmeared();
    esmearing = p->EnergySmearing();
    eerrsmearing = p->EnergyErrSmearing();
    epresmearscale = p->EnergyPreSmearScale();
    escale = p->EnergyScale();    
    idmva = p->IdMva();
    hcalisodr03 = p->HcalTowerSumEtDr03();
    ecalisodr03 = p->EcalRecHitIsoDr03();        
    trkisohollowdr03 = p->HollowConeTrkIsoDr03();
    hcalisodr04 = p->HcalTowerSumEtDr04();
    ecalisodr04 = p->EcalRecHitIsoDr04();        
    trkisohollowdr04 = p->HollowConeTrkIsoDr04();
    
    passeleveto = PhotonTools::PassElectronVetoConvRecovery(p, els, convs, bs);  
    
    const Vertex *vtx = vtxCol->At(0);
    if (p->HasPV()) vtx = p->PV();    
    if ( realVtx ) vtx = realVtx;

    UInt_t wVtxInd = 0;
    
    trackiso1 = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0,
						  trackCol, NULL, NULL,
						  (!applyElectronVeto ? els : NULL) );
    //Question Ming:whyfPV->At(0) instead of selected vertex using ranking method?
    
    // track iso worst vtx
    trackiso2 = IsolationTools::CiCTrackIsolation(p,vtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0,
						  trackCol, &wVtxInd,vtxCol,
						  (!applyElectronVeto ? els : NULL) );
    combiso1 = ecalisodr03+hcalisodr04+trackiso1 - 0.17*rho;
    combiso2 = ecalisodr04+hcalisodr04+trackiso2 - 0.52*rho;  
    
    
    // -----------------------------------------------------
    // PF-CiC4 Debug Stuff
    std::vector<double> debugVals;
    PhotonTools::PassCiCPFIsoSelection(p, vtx, fPFCands, vtxCol, rho, 20., dor9rescale,p0b,p1b,p0e,p1e,&debugVals);
    if( debugVals.size() == 13 ) {
      pfcic4_tIso1   = debugVals[0];
      pfcic4_tIso2   = debugVals[1];
      pfcic4_tIso3   = debugVals[2];
      pfcic4_covIEtaIEta   = debugVals[3];
      pfcic4_HoE   = debugVals[4];
      pfcic4_R9   = debugVals[5];
      pfcic4_wVtxInd   = debugVals[6];
      pfcic4_ecalIso3   = debugVals[7];
      pfcic4_ecalIso4   = debugVals[8];
      pfcic4_trackIsoSel03   = debugVals[9];
      pfcic4_trackIsoWorst04   = debugVals[10];
      pfcic4_combIso1   = debugVals[11];
      pfcic4_combIso2   = debugVals[12];
    }
    // -----------------------------------------------------
    //id mva
    //2011
    idmva_tIso1abs=combiso1;
    idmva_tIso2abs=combiso2;
    idmva_tIso3abs=trackiso1;
    idmva_absIsoEcal=ecalisodr03;
    idmva_absIsoHcal=hcalisodr04;
    //2012  
    idmva_CoviEtaiPhi=p->SCluster()->Seed()->CoviEtaiPhi();
    idmva_s4ratio=p->S4Ratio();
    idmva_GammaIso=IsolationTools::PFGammaIsolation(p,0.3,0,fPFCands);
    idmva_ChargedIso_selvtx=IsolationTools::PFChargedIsolation(p,vtx,0.3,0.,fPFCands);
    idmva_ChargedIso_0p2_selvtx=IsolationTools::PFChargedIsolation(p,vtx,0.2,0.,fPFCands);
    idmva_ChargedIso_worstvtx=IsolationTools::PFChargedIsolation(p,vtx,0.3,0.,fPFCands,&wVtxInd,vtxCol);
    idmva_PsEffWidthSigmaRR=p->EffSigmaRR();
  }
  else {
    hasphoton = kFALSE;
    index = (UInt_t) -1;
    e = -99.;
    pt = -99.;
    eta = -99.;
    phi = -99.;
    r9 = b->E3x3()/s->RawEnergy();
    e3x3 = b->E3x3();
    e5x5 = b->E5x5();
    hovere = ele->HadronicOverEm();
    sigietaieta = ele->CoviEtaiEta();      
    phcat = -99;    
    eerr = -99.;
    eerrsmeared = -99.;
    esmearing = 0.;
    idmva = -99.;
    hcalisodr03 = -99;
    ecalisodr03 = -99;
    trkisohollowdr03 = -99;
    hcalisodr04 = -99;
    ecalisodr04 = -99;
    trkisohollowdr04 = -99;        
    trackiso1 = -99.;
    trackiso2 = -99.;
    combiso1 = -99.;
    combiso2 = -99.;
  }

  if(dor9rescale &&  !isdata){
    if(s->AbsEta()<1.5){
      r9 = p0b + p1b * r9;
    }else{
      r9 = p0e + p1e * r9;
    }
  }
  
  // TODO: fix the bug with supercluster index
  scindex = PhotonTreeWriter::IndexOfNearestSuperClusterInCollection(s, scs);
  /// DEBUG
  // cout << "JV: s, scs, index: " << s << ", " << scs << ", " << scindex << endl;
  sce = s->Energy();
  scrawe = s->RawEnergy();
  scpse = s->PreshowerEnergy();
  scpssigmaxx = s->PsEffWidthSigmaXX();
  scpssigmayy = s->PsEffWidthSigmaYY();
  sceta = s->Eta();
  scphi = s->Phi();
  scnclusters = s->ClusterSize();
  scnhits = s->NHits();
  scetawidth = -99.;
  scphiwidth = -99.;
  if (p) {
    scetawidth = p->EtaWidth();
    scphiwidth = p->PhiWidth();
  }
  else {
    scetawidth = s->EtaWidth();
    scphiwidth = s->PhiWidth();        
  }
  isbarrel = (s->AbsEta()<1.5);
  isr9reco = (isbarrel && r9>0.94) || (!isbarrel && r9>0.95);
  isr9cat = (r9>0.94);
  
  eseed = b->Energy();      
  etaseed = b->Eta();
  phiseed = b->Phi();
  ietaseed = b->IEta();
  iphiseed = b->IPhi();
  ixseed = b->IX();
  iyseed = b->IY();
  etacryseed = b->EtaCry();
  phicryseed = b->PhiCry();
  xcryseed = b->XCry();
  ycryseed = b->YCry();
  thetaaxisseed = b->ThetaAxis();
  phiaxisseed = b->PhiAxis();
  sigietaietaseed = TMath::Sqrt(b->CoviEtaiEta());
  sigiphiphiseed = TMath::Sqrt(b->CoviPhiiPhi());
  if (isnan(sigiphiphiseed)) sigiphiphiseed = -99.;
  covietaiphiseed = b->CoviEtaiPhi();
  if (isnan(covietaiphiseed)) covietaiphiseed = -99.;
  e3x3seed = b->E3x3();
  e5x5seed = b->E5x5();
  emaxseed = b->EMax();
  e2ndseed = b->E2nd();
  etopseed = b->ETop();
  ebottomseed = b->EBottom();
  eleftseed = b->ELeft();
  erightseed = b->ERight();
  e1x3seed = b->E1x3();
  e3x1seed = b->E3x1();
  e1x5seed = b->E1x5();
  e2x2seed = b->E2x2();
  e4x4seed = b->E4x4();
  e2x5maxseed = b->E2x5Max();
  e2x5topseed = b->E2x5Top();
  e2x5bottomseed = b->E2x5Bottom();
  e2x5leftseed = b->E2x5Left();
  e2x5rightseed = b->E2x5Right();
  xseedseed = b->Pos().X();
  yseedseed = b->Pos().Y();
  zseedseed = b->Pos().Z();
  nhitsseed = b->NHits();

  if (b2) {
    ebc2 = b2->Energy();      
    etabc2 = b2->Eta();
    phibc2 = b2->Phi();
    ietabc2 = b2->IEta();
    iphibc2 = b2->IPhi();
    ixbc2 = b2->IX();
    iybc2 = b2->IY();
    etacrybc2 = b2->EtaCry();
    phicrybc2 = b2->PhiCry();
    xcrybc2 = b2->XCry();
    ycrybc2 = b2->YCry();
    thetaaxisbc2 = b2->ThetaAxis();
    phiaxisbc2 = b2->PhiAxis();
    sigietaietabc2 = TMath::Sqrt(b2->CoviEtaiEta());
    sigiphiphibc2 = TMath::Sqrt(b2->CoviPhiiPhi());
    if (isnan(sigiphiphibc2)) sigiphiphibc2 = -99.;
    covietaiphibc2 = b2->CoviEtaiPhi();
    if (isnan(covietaiphibc2)) covietaiphibc2 = -99.;
    e3x3bc2 = b2->E3x3();
    e5x5bc2 = b2->E5x5();
    emaxbc2 = b2->EMax();
    e2ndbc2 = b2->E2nd();
    etopbc2 = b2->ETop();
    ebottombc2 = b2->EBottom();
    eleftbc2 = b2->ELeft();
    erightbc2 = b2->ERight();
    e1x3bc2 = b2->E1x3();
    e3x1bc2 = b2->E3x1();
    e1x5bc2 = b2->E1x5();
    e2x2bc2 = b2->E2x2();
    e4x4bc2 = b2->E4x4();
    e2x5maxbc2 = b2->E2x5Max();
    e2x5topbc2 = b2->E2x5Top();
    e2x5bottombc2 = b2->E2x5Bottom();
    e2x5leftbc2 = b2->E2x5Left();
    e2x5rightbc2 = b2->E2x5Right();
    xbc2bc2 = b2->Pos().X();
    ybc2bc2 = b2->Pos().Y();
    zbc2bc2 = b2->Pos().Z();      
    nhitsbc2= b2->NHits();
  }
  else {
    ebc2 = 0;
    etabc2 = 0;
    phibc2 = 0;
    ietabc2 = 0;
    iphibc2 = 0;
    ixbc2 = 0;
    iybc2 = 0;
    etacrybc2 = 0;
    phicrybc2 = 0;
    xcrybc2 = 0;
    ycrybc2 = 0;
    thetaaxisbc2 = 0;
    phiaxisbc2 = 0;
    sigietaietabc2 = 0;
    sigiphiphibc2 = 0;
    covietaiphibc2 = 0;
    e3x3bc2 = 0;
    e5x5bc2 = 0;
    emaxbc2 = 0;
    e2ndbc2 = 0;
    etopbc2 = 0;
    ebottombc2 = 0;
    eleftbc2 = 0;
    erightbc2 = 0;
    e1x3bc2 = 0;
    e3x1bc2 = 0;
    e1x5bc2 = 0;
    e2x2bc2 = 0;
    e4x4bc2 = 0;
    e2x5maxbc2 = 0;
    e2x5topbc2 = 0;
    e2x5bottombc2 = 0;
    e2x5leftbc2 = 0;
    e2x5rightbc2 = 0;
    xbc2bc2 = 0;
    ybc2bc2 = 0;
    zbc2bc2 = 0;
    nhitsbc2 = 0;
  }

  if (bclast) {
    ebclast = bclast->Energy();      
    etabclast = bclast->Eta();
    phibclast = bclast->Phi();
    ietabclast = bclast->IEta();
    iphibclast = bclast->IPhi();
    ixbclast = bclast->IX();
    iybclast = bclast->IY();
    etacrybclast = bclast->EtaCry();
    phicrybclast = bclast->PhiCry();
    xcrybclast = bclast->XCry();
    ycrybclast = bclast->YCry();
    thetaaxisbclast = bclast->ThetaAxis();
    phiaxisbclast = bclast->PhiAxis();
    sigietaietabclast = TMath::Sqrt(bclast->CoviEtaiEta());
    sigiphiphibclast = TMath::Sqrt(bclast->CoviPhiiPhi());
    if (isnan(sigiphiphibclast)) sigiphiphibclast = -99.;
    covietaiphibclast = bclast->CoviEtaiPhi();
    if (isnan(covietaiphibclast)) covietaiphibclast = -99.;
    e3x3bclast = bclast->E3x3();
    nhitsbclast = bclast->NHits();
  }
  else {
    ebclast = 0;
    etabclast = 0;
    phibclast = 0;
    ietabclast = 0;
    iphibclast = 0;
    ixbclast = 0;
    iybclast = 0;
    etacrybclast = 0;
    phicrybclast = 0;
    xcrybclast = 0;
    ycrybclast = 0;
    thetaaxisbclast = 0;
    phiaxisbclast = 0;
    sigietaietabclast = 0;
    sigiphiphibclast = 0;
    covietaiphibclast = 0;
    e3x3bclast = 0;
    nhitsbclast = 0;
  }

  if (bclast2) {
    ebclast2 = bclast2->Energy();      
    etabclast2 = bclast2->Eta();
    phibclast2 = bclast2->Phi();
    ietabclast2 = bclast2->IEta();
    iphibclast2 = bclast2->IPhi();
    ixbclast2 = bclast2->IX();
    iybclast2 = bclast2->IY();
    etacrybclast2 = bclast2->EtaCry();
    phicrybclast2 = bclast2->PhiCry();
    xcrybclast2 = bclast2->XCry();
    ycrybclast2 = bclast2->YCry();
    thetaaxisbclast2 = bclast2->ThetaAxis();
    phiaxisbclast2 = bclast2->PhiAxis();
    sigietaietabclast2 = TMath::Sqrt(bclast2->CoviEtaiEta());
    sigiphiphibclast2 = TMath::Sqrt(bclast2->CoviPhiiPhi());
    if (isnan(sigiphiphibclast2)) sigiphiphibclast2 = -99.;
    covietaiphibclast2 = bclast2->CoviEtaiPhi();
    if (isnan(covietaiphibclast2)) covietaiphibclast2 = -99.;
    e3x3bclast2 = bclast2->E3x3();
    nhitsbclast2 = bclast2->NHits();
  }
  else {
    ebclast2 = 0;
    etabclast2 = 0;
    phibclast2 = 0;
    ietabclast2 = 0;
    iphibclast2 = 0;
    ixbclast2 = 0;
    iybclast2 = 0;
    etacrybclast2 = 0;
    phicrybclast2 = 0;
    xcrybclast2 = 0;
    ycrybclast2 = 0;
    thetaaxisbclast2 = 0;
    phiaxisbclast2 = 0;
    sigietaietabclast2 = 0;
    sigiphiphibclast2 = 0;
    covietaiphibclast2 = 0;
    e3x3bclast2 = 0;
    nhitsbclast2 = 0;
  }

  for (UInt_t iclus=0; iclus<NClus; ++iclus) {
    if (fillclusterarrays && iclus < s->ClusterSize() ) {
      const BasicCluster *ib =s->Cluster(iclus);
      
      ebcs[iclus] = ib->Energy();      
      etabcs[iclus] = ib->Eta();
      phibcs[iclus] = ib->Phi();
      ietabcs[iclus] = ib->IEta();
      iphibcs[iclus] = ib->IPhi();
      ixbcs[iclus] = ib->IX();
      iybcs[iclus] = ib->IY();
      etacrybcs[iclus] = ib->EtaCry();
      phicrybcs[iclus] = ib->PhiCry();
      xcrybcs[iclus] = ib->XCry();
      ycrybcs[iclus] = ib->YCry();
      sigietaietabcs[iclus] = TMath::Sqrt(ib->CoviEtaiEta());
      sigiphiphibcs[iclus] = TMath::Sqrt(ib->CoviPhiiPhi());
      covietaiphibcs[iclus] = ib->CoviEtaiPhi();
      sigetaetabcs[iclus] = TMath::Sqrt(ib->CovEtaEta());
      sigphiphibcs[iclus] = TMath::Sqrt(ib->CovPhiPhi());
      covetaphibcs[iclus] = ib->CovEtaPhi();
      e3x3bcs[iclus] = ib->E3x3();
      e5x5bcs[iclus] = ib->E5x5();
      emaxbcs[iclus] = ib->EMax();
      e2ndbcs[iclus] = ib->E2nd();
      etopbcs[iclus] = ib->ETop();
      ebottombcs[iclus] = ib->EBottom();
      eleftbcs[iclus] = ib->ELeft();
      erightbcs[iclus] = ib->ERight();
      e1x3bcs[iclus] = ib->E1x3();
      e3x1bcs[iclus] = ib->E3x1();
      e1x5bcs[iclus] = ib->E1x5();
      e2x2bcs[iclus] = ib->E2x2();
      e4x4bcs[iclus] = ib->E4x4();
      e2x5maxbcs[iclus] = ib->E2x5Max();
      e2x5topbcs[iclus] = ib->E2x5Top();
      e2x5bottombcs[iclus] = ib->E2x5Bottom();
      e2x5leftbcs[iclus] = ib->E2x5Left();
      e2x5rightbcs[iclus] = ib->E2x5Right();  
      nhitsbcs[iclus]= ib->NHits();
    }
    else {
      ebcs[iclus] = -999;
      etabcs[iclus] = -999;
      phibcs[iclus] = -999;
      ietabcs[iclus] = -999;
      iphibcs[iclus] = -999;
      ixbcs[iclus] = -999;
      iybcs[iclus] = -999;
      etacrybcs[iclus] = -999;
      phicrybcs[iclus] = -999;
      xcrybcs[iclus] = -999;
      ycrybcs[iclus] = -999;
      sigietaietabcs[iclus] = -999;
      sigiphiphibcs[iclus] = -999;
      covietaiphibcs[iclus] = -999;
      sigetaetabcs[iclus] = -999;
      sigphiphibcs[iclus] = -999;
      covetaphibcs[iclus] = -999;          
      e3x3bcs[iclus] = -999;
      e5x5bcs[iclus] = -999;
      emaxbcs[iclus] = -999;
      e2ndbcs[iclus] = -999;
      etopbcs[iclus] = -999;
      ebottombcs[iclus] = -999;
      eleftbcs[iclus] = -999;
      erightbcs[iclus] = -999;
      e1x3bcs[iclus] = -999;
      e3x1bcs[iclus] = -999;
      e1x5bcs[iclus] = -999;
      e2x2bcs[iclus] = -999;
      e4x4bcs[iclus] = -999;
      e2x5maxbcs[iclus] = -999;
      e2x5topbcs[iclus] = -999;
      e2x5bottombcs[iclus] = -999;
      e2x5leftbcs[iclus] = -999;
      e2x5rightbcs[iclus] = -999;
      nhitsbcs[iclus] = -999;
    }
  }

  for (UInt_t iclus=0; iclus<NClus; ++iclus) {
    if (fillclusterarrays && pfsc && iclus < pfsc->ClusterSize() ) {
      const BasicCluster *ib =pfsc->Cluster(iclus);
      
      epfbcs[iclus] = ib->Energy();      
      etapfbcs[iclus] = ib->Eta();
      phipfbcs[iclus] = ib->Phi();
      ietapfbcs[iclus] = ib->IEta();
      iphipfbcs[iclus] = ib->IPhi();
      ixpfbcs[iclus] = ib->IX();
      iypfbcs[iclus] = ib->IY();
      etacrypfbcs[iclus] = ib->EtaCry();
      phicrypfbcs[iclus] = ib->PhiCry();
      xcrypfbcs[iclus] = ib->XCry();
      ycrypfbcs[iclus] = ib->YCry();
      sigietaietapfbcs[iclus] = TMath::Sqrt(ib->CoviEtaiEta());
      sigiphiphipfbcs[iclus] = TMath::Sqrt(ib->CoviPhiiPhi());
      covietaiphipfbcs[iclus] = ib->CoviEtaiPhi();
      sigetaetapfbcs[iclus] = TMath::Sqrt(ib->CovEtaEta());
      sigphiphipfbcs[iclus] = TMath::Sqrt(ib->CovPhiPhi());
      covetaphipfbcs[iclus] = ib->CovEtaPhi();
      e3x3pfbcs[iclus] = ib->E3x3();
      e5x5pfbcs[iclus] = ib->E5x5();
      emaxpfbcs[iclus] = ib->EMax();
      e2ndpfbcs[iclus] = ib->E2nd();
      etoppfbcs[iclus] = ib->ETop();
      ebottompfbcs[iclus] = ib->EBottom();
      eleftpfbcs[iclus] = ib->ELeft();
      erightpfbcs[iclus] = ib->ERight();
      e1x3pfbcs[iclus] = ib->E1x3();
      e3x1pfbcs[iclus] = ib->E3x1();
      e1x5pfbcs[iclus] = ib->E1x5();
      e2x2pfbcs[iclus] = ib->E2x2();
      e4x4pfbcs[iclus] = ib->E4x4();
      e2x5maxpfbcs[iclus] = ib->E2x5Max();
      e2x5toppfbcs[iclus] = ib->E2x5Top();
      e2x5bottompfbcs[iclus] = ib->E2x5Bottom();
      e2x5leftpfbcs[iclus] = ib->E2x5Left();
      e2x5rightpfbcs[iclus] = ib->E2x5Right();  
      nhitspfbcs[iclus]= ib->NHits();
    }
    else {
      epfbcs[iclus] = -999;
      etapfbcs[iclus] = -999;
      phipfbcs[iclus] = -999;
      ietapfbcs[iclus] = -999;
      iphipfbcs[iclus] = -999;
      ixpfbcs[iclus] = -999;
      iypfbcs[iclus] = -999;
      etacrypfbcs[iclus] = -999;
      phicrypfbcs[iclus] = -999;
      xcrypfbcs[iclus] = -999;
      ycrypfbcs[iclus] = -999;
      sigietaietapfbcs[iclus] = -999;
      sigiphiphipfbcs[iclus] = -999;
      covietaiphipfbcs[iclus] = -999;
      sigetaetapfbcs[iclus] = -999;
      sigphiphipfbcs[iclus] = -999;
      covetaphipfbcs[iclus] = -999;          
      e3x3pfbcs[iclus] = -999;
      e5x5pfbcs[iclus] = -999;
      emaxpfbcs[iclus] = -999;
      e2ndpfbcs[iclus] = -999;
      etoppfbcs[iclus] = -999;
      ebottompfbcs[iclus] = -999;
      eleftpfbcs[iclus] = -999;
      erightpfbcs[iclus] = -999;
      e1x3pfbcs[iclus] = -999;
      e3x1pfbcs[iclus] = -999;
      e1x5pfbcs[iclus] = -999;
      e2x2pfbcs[iclus] = -999;
      e4x4pfbcs[iclus] = -999;
      e2x5maxpfbcs[iclus] = -999;
      e2x5toppfbcs[iclus] = -999;
      e2x5bottompfbcs[iclus] = -999;
      e2x5leftpfbcs[iclus] = -999;
      e2x5rightpfbcs[iclus] = -999;
      nhitspfbcs[iclus] = -999;
    }
  }

  for (UInt_t iclus=0; iclus<100; ++iclus) {
    if (fillclusterarrays && pfsc && iclus < pfsc->NPsClusts() ) {
      const PsCluster *ib = pfsc->PsClust(iclus);
      
      epsc[iclus] = ib->Energy();
      etapsc[iclus] = ib->Eta();
      phipsc[iclus] = ib->Phi();
      planepsc[iclus] = ib->PsPlane();
    }
    else {
      epsc[iclus] = -999;
      etapsc[iclus] = -999;
      phipsc[iclus] = -999;
      planepsc[iclus] = 0;      
    }
  }

  //initialize photon energy corrections if needed
  /*if (!PhotonFix::initialised()) {
    PhotonFix::initialise("4_2",std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFix.dat")).Data()));  
  }*/ 

  phfixph.setup(e,sceta,scphi,r9);
  phfixele.setup(e,sceta,scphi,r9);      
  
  const Float_t dval = -99.;
  ecor = phfixph.fixedEnergy();
  ecorerr = phfixph.sigmaEnergy();
  ecorele = phfixele.fixedEnergy();
  ecoreleerr = phfixele.sigmaEnergy();
  if (phfixph.isbarrel()) {
   etac = phfixph.etaC();
   etas = phfixph.etaS();
   etam = phfixph.etaM();
   phic = phfixph.phiC();
   phis = phfixph.phiS();
   phim = phfixph.phiM();
   xz = dval;
   xc = dval;
   xs = dval;
   xm = dval;
   yz = dval;
   yc = dval;
   ys = dval;
   ym = dval;
  }
  else {
   etac = dval;
   etas = dval;
   etam = dval;
   phic = dval;
   phis = dval;
   phim = dval;
   xz = phfixph.xZ();
   xc = phfixph.xC();
   xs = phfixph.xS();
   xm = phfixph.xM();
   yz = phfixph.yZ();
   yc = phfixph.yC();
   ys = phfixph.yS();
   ym = phfixph.yM();
  }   
  
  if (c) {
    hasconversion = kTRUE;
    convp = c->P();
    convpt = c->Pt();
    conveta = c->Eta();
    convphi = c->Phi();
    ThreeVector dirconvsc = ThreeVector(s->Point()) - c->Position();
    convdeta = c->Eta() - dirconvsc.Eta();
    convdphi = MathUtils::DeltaPhi(c->Phi(),dirconvsc.Phi());
    convvtxrho = c->Position().Rho();
    convvtxz =   c->Position().Z();
    convvtxphi = c->Position().Phi();        
   
    const StableData *leadsd = dynamic_cast<const StableData*>(c->DaughterDat(0));
    const StableData *trailsd = dynamic_cast<const StableData*>(c->DaughterDat(1));
    if (leadsd->Pt()<trailsd->Pt()) {
      const StableData *sdtmp = leadsd;
      leadsd = trailsd;
      trailsd = sdtmp;
    }
    
    const Track *leadtrack = dynamic_cast<const StableParticle*>(leadsd->Original())->Trk();
    const Track *trailtrack = dynamic_cast<const StableParticle*>(trailsd->Original())->Trk();
    
    convleadpt = leadsd->Pt();
    convtrailpt = trailsd->Pt();
    convleadtrackpt = leadtrack->Pt();
    convleadtrackalgo = leadtrack->Algo();
    if (convleadtrackalgo==29) convleadtrackalgos=2; //gsf track
    else if (convleadtrackalgo==15 ||convleadtrackalgo==16) convleadtrackalgos=1; //ecal-seeded track
    else convleadtrackalgos = 0; //std iterative track
    convleadtrackcharge = leadtrack->Charge();
    convtrailtrackpt = trailtrack->Pt();
    convtrailtrackalgo = trailtrack->Algo();
    if (convtrailtrackalgo==29) convtrailtrackalgos=2; //gsf track
    else if (convtrailtrackalgo==15 ||convtrailtrackalgo==16) convtrailtrackalgos=1; //ecal-seeded track
    else convtrailtrackalgos = 0; //std iterative track
    trailtrackcharge = trailtrack->Charge();      
  }
  else {
    hasconversion = kFALSE;
    convp = -99.;
    convpt = -99.;
    conveta = -99.;
    convphi = -99.;
    convdeta = -99.;
    convdphi = -99.;
    convvtxrho = -99.;
    convvtxz = -999.;
    convvtxphi = -99.;        
    convleadpt = -99.;
    convtrailpt = -99.;
    convleadtrackpt = -99.;
    convleadtrackalgo = -99;
    convleadtrackalgos = -99;
    convleadtrackcharge = 0;
    convtrailtrackpt = -99.;
    convtrailtrackalgo = -99;
    convtrailtrackalgos = -99;      
    trailtrackcharge = 0;      
  }
  
  //electron quantities
  if (ele) {
    haselectron = kTRUE;
    eleisecaldriven = ele->IsEcalDriven();
    eleistrackerdriven = ele->IsTrackerDriven();
    elee = ele->E();
    elept = ele->Pt();
    eleeta = ele->Eta();
    elephi = ele->Phi();
    elecharge = ele->Charge();
    elefbrem = ele->FBrem();
    eledeta = ele->DeltaEtaSuperClusterTrackAtVtx();
    eledphi = ele->DeltaPhiSuperClusterTrackAtVtx();
    elep = s->Energy()/ele->ESuperClusterOverP();
    elepin = ele->PIn();
    elepout = ele->POut();        
  }
  else {
    haselectron = kFALSE;
    eleisecaldriven = kFALSE;
    eleistrackerdriven = kFALSE;
    elee = -99.;
    elept = -99.;
    eleeta = -99.;
    elephi = -99.;
    elecharge = -99;
    elefbrem = -99.;
    eledeta = -99.;
    eledphi = -99.;
    elep = -99.;
    elepin = -99.;
    elepout = -99.;
  }
  
  //pf supercluster quantities
  if (pfsc) {
    haspfsc = kTRUE;
    pfsce = pfsc->Energy();
    pfscrawe = pfsc->RawEnergy();
    pfsceta = pfsc->Eta();
    pfscphi = pfsc->Phi();        
    pfscnclusters = pfsc->NClusters();
    pfscnhits = pfsc->NHits();
    pfscetawidth = pfsc->EtaWidth();
    pfscphiwidth = pfsc->PhiWidth();
    pfscnpsclusters = pfsc->NPsClusts();   
  }
  else {
    haspfsc = kFALSE;
    pfsce = -99.;
    pfscrawe = -99.;
    pfsceta = -99.;
    pfscphi = -99.;
    pfscnclusters = 0;
    pfscnhits = 0;
    pfscetawidth = -99.;
    pfscphiwidth = -99.;
    pfscnpsclusters = 0;
  }
  
  genz = -99.;
  if (m) {
    ispromptgen = kTRUE;
    gene = m->E();
    genpt = m->Pt();
    geneta = m->Eta();
    genphi = m->Phi();
    const MCParticle *mm = m->DistinctMother();
    if (mm) genz = mm->DecayVertex().Z();
    pdgid = m->PdgId();
    if (mm) motherpdgid = mm->PdgId();
    else motherpdgid = -99;
  }
  else {
    ispromptgen = kFALSE;
    gene = -99.;
    genpt = -99.;
    geneta = -99.;
    genphi = -99.;
    pdgid = -99;
    motherpdgid = -99;
  }
}

void PhotonTreeWriterVtx::SetVars(const Vertex *v, const Photon *p1, const Photon *p2, const PFCandidateCol *pfcands, Int_t idx, Int_t numvtx, Float_t genvtxz)
{
 
  //printf("start\n");
  
  n = idx;
  nvtx = numvtx;
  zgen = genvtxz;
  
  x = v->X();
  y = v->Y();
  z = v->Z();
  
  Double_t dsumpt = 0.;
  Double_t dsumptsq = 0.;
  
  nchalltoward = 0;
  nchalltransverse = 0;
  nchallaway = 0;
  nchcuttoward = 0;
  nchcuttransverse = 0;
  nchcutaway = 0;  
  
  ThreeVector vtxmom;
  
  //printf("mom\n");
  FourVectorM diphoton = p1->MomVtx(v->Position()) + p2->MomVtx(v->Position());
  //printf("done mom\n");
  ptgg = diphoton.Pt();
  phigg = diphoton.Phi();
  etagg = diphoton.Eta();
  mgg = diphoton.M();
  pxgg = diphoton.Px();
  pygg = diphoton.Py();
  pzgg = diphoton.Pz();
  
  //printf("loop\n");
  
  for (UInt_t i = 0; i<pfcands->GetEntries(); ++i) {
    const PFCandidate *pfc = pfcands->At(i);
    if (pfc->PFType()!=PFCandidate::eHadron || !pfc->HasTrackerTrk()) continue;
    if (TMath::Abs( pfc->TrackerTrk()->DzCorrected(*v) ) > 0.2) continue;
    if (TMath::Abs( pfc->TrackerTrk()->D0Corrected(*v) ) > 0.1) continue;
    
    vtxmom += ThreeVector(pfc->Px(),pfc->Py(),pfc->Pz());
    
    dsumpt += pfc->Pt();
    dsumptsq += pfc->Pt()*pfc->Pt();
    
    Double_t dphi = TMath::Abs(MathUtils::DeltaPhi(*pfc,diphoton));
    if (dphi<(TMath::Pi()/3.0)) {
      ++nchalltoward;
      if (pfc->Pt()>0.5) ++nchcuttoward;
    }
    else if (dphi>(2.0*TMath::Pi()/3.0)) {
      ++nchallaway;
      if (pfc->Pt()>0.5) ++nchcutaway;      
    }
    else {
      ++nchalltransverse;
      if (pfc->Pt()>0.5) ++nchcuttransverse;            
    }
     
  }
  
  //printf("doneloop\n");

    
  sumpt = dsumpt;
  sumptsq = dsumptsq;
  
  pt = vtxmom.Rho();
  phi = vtxmom.Phi();
  eta = vtxmom.Eta();
  px = vtxmom.X();
  py = vtxmom.Y();
  pz = vtxmom.Z();
  
  //printf("done\n");
  
  return;
  
}


//_____________________________________________________________________________
void PhotonTreeWriter::ApplyLeptonTag(const Photon *phHard,
                                      const Photon *phSoft,
                                      const Vertex *selvtx)
{
  // perform flavor-based lepton tagging (used before the legacy paper of 2013)
  // the diphoton event record will have one more entry; i.e. leptonTag
  // leptonTag = -1   -> lepton-taggng was swicthed off
  //           =  0   -> event tagged as 'non-lepton-event'
  //           = +1   -> event tagged as muon-event
  //           = +2   -> event tagged as electron-event
  fDiphotonEvent->leptonTag = 0;
  Int_t closestVtx = 0;
  const Muon *muon = GetLeptonTagMuon(phHard, phSoft);
  if (muon) { 
    // muon tagged
    fDiphotonEvent->leptonTag = 2;
    SetLeptonTagMuonVars(phHard, phSoft, muon);
  } 
  
  const Electron *electron = 0;
  if (fDiphotonEvent->leptonTag < 1) {
    electron = GetLeptonTagElectron(phHard, phSoft);
  }
  
  if (electron && 
      MassOfPairIsWithinWindowAroundMZ(phHard, electron, 10) == false &&
      MassOfPairIsWithinWindowAroundMZ(phSoft, electron, 10) == false) {
    // electron tagged
    fDiphotonEvent->leptonTag = 1;
    
    fDiphotonEvent-> elePt = electron->Pt();
    fDiphotonEvent-> eleEta = electron->Eta();
    fDiphotonEvent-> elePhi = electron->Phi();
    fDiphotonEvent-> eleSCEta = electron->SCluster()->Eta();
    fDiphotonEvent-> eleIso1 = (electron->TrackIsolationDr03() + electron->EcalRecHitIsoDr03() + electron->HcalTowerSumEtDr03() - fPileUpDen->At(0)->RhoRandomLowEta() * TMath::Pi() * 0.3 * 0.3)/fDiphotonEvent-> elePt;
    
    fDiphotonEvent-> eleIso2 = -99.;
    
    if ( fDoSynching ) {
      Double_t distVtx = 999.0;
      for(UInt_t nv=0; nv<fPV->GetEntries(); nv++){
        double dz = TMath::Abs(electron->GsfTrk()->DzCorrected(*fPV->At(nv)));
        if(dz < distVtx) {
          distVtx    = dz;
          closestVtx = nv;
        }
      }
      fDiphotonEvent-> eleIdMva = fElectronIDMVA->MVAValue(electron, fPV->At(closestVtx));
    }
    
    
    fDiphotonEvent-> eleIso3 = (electron->TrackIsolationDr03() + electron->EcalRecHitIsoDr03() + electron->HcalTowerSumEtDr03() - fPileUpDen->At(0)->RhoLowEta() * TMath::Pi() * 0.3 * 0.3)/fDiphotonEvent-> elePt;
    fDiphotonEvent-> eleIso4 = (electron->TrackIsolationDr03() + electron->EcalRecHitIsoDr03() + electron->HcalTowerSumEtDr03() - fPileUpDen->At(0)->Rho() * TMath::Pi() * 0.3 * 0.3)/fDiphotonEvent-> elePt;
    fDiphotonEvent-> eleDist = electron->ConvPartnerDist();
    fDiphotonEvent-> eleDcot = electron->ConvPartnerDCotTheta();
    fDiphotonEvent-> eleCoviee = electron->CoviEtaiEta();
    fDiphotonEvent-> eleDphiin = TMath::Abs(electron->DeltaPhiSuperClusterTrackAtVtx());
    fDiphotonEvent-> eleDetain = TMath::Abs(electron->DeltaEtaSuperClusterTrackAtVtx());
    fDiphotonEvent-> eleDR1 = MathUtils::DeltaR(electron,phHard);
    fDiphotonEvent-> eleDR2 = MathUtils::DeltaR(electron,phSoft);
    fDiphotonEvent-> eleMass1 = (phHard->Mom()+electron->Mom()).M();
    fDiphotonEvent-> eleMass2 = (phSoft->Mom()+electron->Mom()).M();
    fDiphotonEvent-> eleNinnerHits =      electron->Trk()->NExpectedHitsInner();
  } // electron tagged
  
  if(false){
    if(fDiphotonEvent->evt==79737729 || fDiphotonEvent->evt== 871378986  || fDiphotonEvent->evt==528937923 || fDiphotonEvent->evt== 261543921){
      printf("ming sync check ele:  run:%d  evt:%d  lumi:%d  leptonTag:%d  numelectrons:%d  idmva:%f  mass:%f\n  elePt:%f  eleEta:%f  eleSCEta:%f  vtx:%d\n",fDiphotonEvent->run,fDiphotonEvent->evt,fDiphotonEvent->lumi,fDiphotonEvent->leptonTag,fLeptonTagElectrons->GetEntries(),fDiphotonEvent->eleIdMva,fDiphotonEvent->mass,fDiphotonEvent->elePt,fDiphotonEvent->eleEta,fDiphotonEvent->eleSCEta,closestVtx);
      //return;
    }
    if(fDiphotonEvent->evt==333643114 || fDiphotonEvent->evt==89022540 || fDiphotonEvent->evt==8983064 || fDiphotonEvent->evt==876316897 || fDiphotonEvent->evt==541603559  || fDiphotonEvent->evt==223740859) {
      printf("ming sync check muon:  run:%d  evt:%d  lumi:%d  leptonTag:%d  numMuons:%d  mass:%f\n  muonPt:%f  muonEta:%f\n\n",fDiphotonEvent->run,fDiphotonEvent->evt,fDiphotonEvent->lumi,fDiphotonEvent->leptonTag,fLeptonTagMuons->GetEntries(),fDiphotonEvent->mass,fDiphotonEvent->muonPt,fDiphotonEvent->muonEta);
      //return;
    }
  }
} // void PhotonTreeWriter::ApplyLeptonTag(..)


//_____________________________________________________________________________
const Muon* PhotonTreeWriter::GetLeptonTagMuon(const Photon *phHard, 
                                               const Photon *phSoft)
{
  // need to have dR > 1 for with respect to both photons ***changed to 0.7 for 2012
  if (fVerbosityLevel > 90) {
    cout << "++ Entering GetLeptonTagMuon ..." << endl
         << "++ fLeptonTagMuons->GetEntries():                    " 
         << fLeptonTagMuons->GetEntries()               << endl;
    const Muon *muon = 0;
    if (fLeptonTagMuons->GetEntries() > 0) {
      muon = fLeptonTagMuons->At(0);
    }
    if (muon) {
      cout << "++ fLeptonTagMuons->At(0):                           " 
           << fLeptonTagMuons->At(0)                      << endl
           << "++ MathUtils::DeltaR(fLeptonTagMuons->At(0), phHard): " 
           << MathUtils::DeltaR(fLeptonTagMuons->At(0), phHard) << endl
           << "++ MathUtils::DeltaR(fLeptonTagMuons->At(0), phSoft): " 
           << MathUtils::DeltaR(fLeptonTagMuons->At(0), phSoft) << endl
           << "++ muon pt :" << muon->Pt () << endl
           << "++ muon eta:" << muon->Eta() << endl
           << "++ muon phi:" << muon->Phi() << endl;
    }
    cout << "++ Exiting GetLeptonTagMuon ..." << endl;
  }
  
  for (unsigned int imu=0; imu<fLeptonTagMuons->GetEntries(); ++imu) {
    const Muon *mu = fLeptonTagMuons->At(imu);
    if (fLeptonTagMuons->GetEntries() > 0                       &&
        fLeptonTagMuons->At(0) != 0                             &&
        MathUtils::DeltaR(mu, phHard) >= 1.0 && 
        MathUtils::DeltaR(mu, phSoft) >= 1.0) {
      return mu;
    }
  }
    
  return 0;
  
} // const Muon* PhotonTreeWriter::GetLeptonTagMuon(..)
  
  
//_____________________________________________________________________________
void PhotonTreeWriter::SetLeptonTagMuonVars(const Photon *phHard, 
                                            const Photon *phSoft,
                                            const Muon *muon)
{
  fDiphotonEvent-> muonPt  = muon->Pt();
  fDiphotonEvent-> muonEta = muon->Eta();
  fDiphotonEvent-> muonPhi = muon->Phi();
  fDiphotonEvent-> muDR1   = MathUtils::DeltaR(muon, phHard);
  fDiphotonEvent-> muDR2   = MathUtils::DeltaR(muon, phSoft);
  
  Float_t combinedIso = (muon->IsoR03SumPt() + 
                         muon->IsoR03EmEt() + 
                         muon->IsoR03HadEt());
  
  Float_t coneArea = TMath::Pi() * 0.3 * 0.3;
  
  Float_t rho1 = fPileUpDen->At(0)->RhoRandomLowEta();
  Float_t rho2 = fPileUpDen->At(0)->RhoRandom();
  Float_t rho3 = fPileUpDen->At(0)->RhoLowEta();
  Float_t rho4 = fPileUpDen->At(0)->Rho();
  
  fDiphotonEvent-> muIso1 = (combinedIso - rho1 * coneArea) / muon->Pt(); 
  fDiphotonEvent-> muIso2 = (combinedIso - rho2 * coneArea) / muon->Pt(); 
  fDiphotonEvent-> muIso3 = (combinedIso - rho3 * coneArea) / muon->Pt(); 
  fDiphotonEvent-> muIso4 = (combinedIso - rho4 * coneArea) / muon->Pt(); 
  
  fDiphotonEvent-> muD0  = TMath::Abs(muon->BestTrk()->D0Corrected(*fPV->At(0)));
  fDiphotonEvent-> muDZ  = TMath::Abs(muon->BestTrk()->DzCorrected(*fPV->At(0)));
  fDiphotonEvent-> muChi2  = muon->GlobalTrk()->Chi2()/muon->GlobalTrk()->Ndof();
  
  fDiphotonEvent-> muNhits = muon->BestTrk()->NHits();
  fDiphotonEvent-> muNpixhits = muon->BestTrk()->NPixelHits();
  fDiphotonEvent-> muNegs = muon->NSegments();
  fDiphotonEvent-> muNMatch = muon->NMatches();
} // void PhotonTreeWriter::SetLeptonTagMuonVars(..)


//_____________________________________________________________________________
const Electron* PhotonTreeWriter::GetLeptonTagElectron(const Photon *phHard, 
                                                       const Photon *phSoft)
{
  // need to have dR > 1 for with respect to both photons ***changed to 0.7 for 2012
  if (fVerbosityLevel > 90) {
    cout << "++ Entering GetLeptonTagElectron ..." << endl
         << "++ fLeptonTagElectrons->GetEntries():                         " 
         <<  fLeptonTagElectrons->GetEntries()                          << endl
         << "++ PhotonTools::ElectronVetoCiC(phHard, fLeptonTagElectrons): " 
         <<  PhotonTools::ElectronVetoCiC(phHard, fLeptonTagElectrons)  << endl
         << "++ PhotonTools::ElectronVetoCiC(phSoft, fLeptonTagElectrons): " 
         <<  PhotonTools::ElectronVetoCiC(phSoft, fLeptonTagElectrons)  << endl
         << "++ PhotonTools::ElectronVetoCiC(phHard, fElectrons):          " 
         <<  PhotonTools::ElectronVetoCiC(phHard, fElectrons)           << endl
         << "++ PhotonTools::ElectronVetoCiC(phSoft, fElectrons):          " 
         <<  PhotonTools::ElectronVetoCiC(phSoft, fElectrons)           << endl;
     if (fLeptonTagElectrons->GetEntries() > 0 && fLeptonTagElectrons->At(0) != 0) {
       cout  << "++ fLeptonTagElectrons->At(0):                                " 
             <<  fLeptonTagElectrons->At(0)                                 << endl
             << "++ MathUtils::DeltaR(fLeptonTagElectrons->At(0), phHard):     " 
             <<  MathUtils::DeltaR(fLeptonTagElectrons->At(0), phHard)  << endl
             << "++ MathUtils::DeltaR(fLeptonTagElectrons->At(0), phSoft):     " 
             <<  MathUtils::DeltaR(fLeptonTagElectrons->At(0), phSoft)  << endl
             << "++ electron pt : " << fLeptonTagElectrons->At(0)->Pt () << endl
             << "++ electron eta: " << fLeptonTagElectrons->At(0)->Eta() << endl
             << "++ electron phi: " << fLeptonTagElectrons->At(0)->Phi() << endl;
     }
     cout << "++ Exiting GetLeptonTagElectron ..." << endl;
  }
  
  //loop over all electrons
  for (unsigned int iele=0; iele<fLeptonTagElectrons->GetEntries(); ++iele) {
    const Electron *ele = fLeptonTagElectrons->At(iele);
    if (fLeptonTagElectrons->GetEntries() > 0                          &&
        fLeptonTagElectrons->At(0) != 0                                &&
        PhotonTools::ElectronVetoCiC(phHard, fElectrons) >= 1.         && 
        PhotonTools::ElectronVetoCiC(phSoft, fElectrons) >= 1.         &&     
        MathUtils::DeltaR(ele, phHard) >= 1.     &&
        MathUtils::DeltaR(ele, phSoft) >= 1.){
      return ele;
    } 
  }
    
  return 0;
  
} // const Electron* PhotonTreeWriter::GetLeptonTagElectron(..)

  
//_____________________________________________________________________________
bool PhotonTreeWriter::MassOfPairIsWithinWindowAroundMZ(
  const Particle * particle1,
  const Particle * particle2,
  Float_t halfWindowSize,
  Float_t MZ
) {
  Float_t mass = (particle1->Mom() + particle2->Mom()).M();
  return TMath::Abs(mass - MZ) < halfWindowSize;
} // bool PhotonTreeWriter::MassOfPairIsWithinWindowAroundMZ(..)


//_____________________________________________________________________________
void PhotonTreeWriter::ApplyVHLepTag(const Photon *phHard,
                                     const Photon *phSoft,
                                     const Vertex *selvtx)
{
  
  // Perform flavor-based lepton tagging (used since the legacy paper of 2013)
  // the diphoton event record will have one more entry; i.e. leptonTag
  // VHLepTag = -1   -> lepton-tagging was swicthed off
  //          =  0   -> event tagged as 'non-lepton event'
  //          = +1   -> event tagged as a low-MET low-S/sqrt(B) event
  //          = +2   -> event tagged as a high-MET/dilepton high-S/sqrt(B) event
  //          = +3   -> event tagged as both a loose and a tight VH leptonic event
  // TODO: Set the selected vertex to the lepton vertex if tagged as a
  //       VH(lep) event.
  
  //printf("ApplyVHLepTag: electrons = %i, softelectrons = %i, muons = %i, softmuons = %i\n",fLeptonTagElectrons->GetEntries(), fLeptonTagSoftElectrons->GetEntries(), fLeptonTagMuons->GetEntries(), fLeptonTagSoftMuons->GetEntries());
  
  if (fVerbosityLevel > 90) {
    cout << "+ Entering ApplyVHLepTag ..." << endl
         << "+ fElectrons->GetEntries(): " << fElectrons->GetEntries() << endl
         << "+ fLeptonTagElectrons->GetEntries(): " 
         << fLeptonTagElectrons->GetEntries() << endl
         << "+ fLeptonTagMuons->GetEntries(): " 
         << fLeptonTagMuons->GetEntries() << endl
         << "+ corrpfmet: " << fDiphotonEvent->corrpfmet << endl;
  }

  fDiphotonEvent->VHLepTag = 0; // non-lepton event
  bool isVHLepLoose = false;
  bool isVHLepTight = false;

  if (VHLepHasDielectron(phHard, phSoft) || 
      VHLepHasDimuon    (phHard, phSoft)) isVHLepTight = true;
  
  if (fDiphotonEvent->leptonTag < 0) {
    ApplyLeptonTag(phHard, phSoft, selvtx);
  }
  
  //printf("check dilepton, tight = %i\n",int(isVHLepTight));

  const Muon *muon = GetLeptonTagMuon(phHard, phSoft);
  if (muon && VHLepNumberOfJets(phHard, phSoft, selvtx, muon) <= 2) {
    if (fDiphotonEvent->corrpfmet > 45.) isVHLepTight = true; // high MET event
    else                                 isVHLepLoose = true; // low MET event
  } // Found a good VH(lep) tag muon.
  
  //printf("check muon, tight  = %i, loose = %i\n",int(isVHLepTight),int(isVHLepLoose));
  
  const Electron *electron = GetLeptonTagElectron(phHard, phSoft);
  if (electron && VHLepNumberOfJets(phHard, phSoft, selvtx, electron) <= 2) {
    if (fDiphotonEvent->corrpfmet > 45.) {
      isVHLepTight = true;      // High MET event.
    } else if (!MassOfPairIsWithinWindowAroundMZ(phHard, electron, 10) &&
               !MassOfPairIsWithinWindowAroundMZ(phSoft, electron, 10)) {
      isVHLepLoose = true;   // low MET event
    } // Low MET event.
    fDiphotonEvent->elePt  = electron->Pt ();
    fDiphotonEvent->eleEta = electron->Eta();
    fDiphotonEvent->elePhi = electron->Phi();
  } // Found a good VH(lep) tag electron.

  //printf("check electron, tight  = %i, loose = %i\n",int(isVHLepTight),int(isVHLepLoose));
  
  if (isVHLepLoose) fDiphotonEvent->VHLepTag += 1;
  if (isVHLepTight) fDiphotonEvent->VHLepTag += 2;

  if (fVerbosityLevel > 90) {
    cout << "+ muon: " <<  muon << endl;
    if (muon) {
      cout << "+ njets(mu): " 
           << VHLepNumberOfJets(phHard, phSoft, selvtx, muon) << endl;
    }         
    cout << "+ electron: " <<  electron << endl;
    if (electron) {
      cout << "+ njets(ele): " 
           << VHLepNumberOfJets(phHard, phSoft, selvtx, electron) << endl
           << "+ MassOfPairIsWithinWindowAroundMZ(phHard, electron, 10)" 
           <<    MassOfPairIsWithinWindowAroundMZ(phHard, electron, 10) << endl         
           << "+ MassOfPairIsWithinWindowAroundMZ(phSoft, electron, 10)" 
           <<    MassOfPairIsWithinWindowAroundMZ(phSoft, electron, 10) << endl;
    }
    cout << "+ isVHLepTight: " << isVHLepTight << endl
         << "+ isVHLepLoose: " << isVHLepLoose << endl
         << "+ VHLepTag: " << fDiphotonEvent->VHLepTag << endl
         << "+ Exiting ApplyVHLepTag ... " << endl;
  } /// fVerbosityLevel > 90
  
} // void PhotonTreeWriter::ApplyVHLepTag(..)


//_____________________________________________________________________________
void PhotonTreeWriter::ApplyVHHadTag(const Photon *phHard,
                                     const Photon *phSoft,
                                     const Vertex *selvtx,
                                     const Jet    *jet1,
                                     const Jet    *jet2)
{
  
  // Perform VH hadronic tagging (used since the legacy paper of 2013)
  // the diphoton event record will have one more entry; i.e. VHHadTag
  // VHHadTag = -1   -> VH(had)-tagging was swicthed off
  //          =  0   -> event tagged as 'non-VH(had) event'
  //          = +1   -> event tagged as a VH(had) event
  // TODO: Should the selected vertex be updated using the vertex of the jets?

  fDiphotonEvent->VHHadTag = 0; // non-VH(had) event
  
  // Calculate |cos(theta*)|
  // Inspired by Globe:
  // https://github.com/h2gglobe/h2gglobe/blob/ae4356ac0d18e6f77da6e0420ab6f5168e353315/PhotonAnalysis/src/PhotonAnalysis.cc#L4369-L4377
  FourVectorM pgg = phHard->Mom() + phSoft->Mom(); // diphoton 4-vector
  FourVectorM pjj = jet1  ->Mom() + jet2  ->Mom(); // dijet 4-vector
  FourVectorM p4   = pjj + pgg;                   // 4-body 4-vector
  TLorentzVector Vstar(p4.X(), p4.Y(), p4.Z(), p4.T());
  TLorentzVector H_Vstar(pgg.X(), pgg.Y(), pgg.Z(), pgg.T());
  H_Vstar.Boost(-Vstar.BoostVector());
  double cosThetaStar = -H_Vstar.CosTheta();
  double absCosThetaStar = TMath::Abs(cosThetaStar);

  if (fDoSynching) {
    fDiphotonEvent->costhetastar = cosThetaStar;
  }

  Float_t &mass       = fDiphotonEvent->mass;
  Float_t reducedMass = mass / 120.;
  UInt_t  nJets       = NumberOfJets(phHard, phSoft, selvtx, 40., 2.4);
  Float_t mjj         = fDiphotonEvent->dijetmass;

  /// See L2007-2013 of the Hgg AN 2013/253 v3
  if (fDiphotonEvent->ptgg > 130. * reducedMass &&
      nJets >= 2 &&
      60. < mjj && mjj < 120. &&
      absCosThetaStar < 0.5) {
    // Tag this event as a VH-hadronic one!
    fDiphotonEvent->VHHadTag = 1;
  }
} // void PhotonTreeWriter::ApplyVHHadTag(..)  


//_____________________________________________________________________________
bool PhotonTreeWriter::VHLepHasDielectron(const Photon *phHard,
                                          const Photon *phSoft) {
  if (fLeptonTagSoftElectrons->GetEntries() < 2) return false;
  
  if (PhotonTools::ElectronVetoCiC(phHard, fElectrons)              < 1. || 
      PhotonTools::ElectronVetoCiC(phSoft, fElectrons)              < 1.) {
    return false;
  }

  vector<UInt_t> goodElectrons;

  // Loop over electrons.
  for (UInt_t iele=0; iele < fLeptonTagSoftElectrons->GetEntries(); ++iele){
    const Electron *ele = fLeptonTagSoftElectrons->At(iele);
    if (MathUtils::DeltaR(ele, phHard) < 1.0) continue;
    if (MathUtils::DeltaR(ele, phSoft) < 1.0) continue;
    goodElectrons.push_back(iele);
  }
  
  if (goodElectrons.size() < 2) return false;

  // Loop over pairs of selected electrons.
  for (UInt_t iiele1 = 0; iiele1 < goodElectrons.size() - 1; ++iiele1) {
    UInt_t iele1 = goodElectrons[iiele1];
    const Electron *ele1 = fLeptonTagSoftElectrons->At(iele1);
    for (UInt_t iiele2 = iiele1 + 1; iiele2 < goodElectrons.size(); ++iiele2) {
      UInt_t iele2 = goodElectrons[iiele2];
      const Electron *ele2 = fLeptonTagSoftElectrons->At(iele2);
      Double_t mass12 = (ele1->Mom() + ele2->Mom()).M();
      if (mass12 < 70. || 110. < mass12) continue;
      if (fVerbosityLevel > 0) {
        cout << "    Found a tagging dielectron!" << endl << flush;
      }
      fDiphotonEvent->ele1Pt  = ele1->Pt ();
      fDiphotonEvent->ele1Eta = ele1->Eta();
      fDiphotonEvent->ele1Phi = ele1->Phi();
      fDiphotonEvent->ele2Pt  = ele2->Pt ();
      fDiphotonEvent->ele2Eta = ele2->Eta();
      fDiphotonEvent->ele2Phi = ele2->Phi();
      return true;
    }
  }
  
  return false;
  
} // bool PhotonTreeWriter::VHLepHasDielectron(..)


//_____________________________________________________________________________
bool PhotonTreeWriter::VHLepHasDimuon(const Photon *phHard,
                                      const Photon *phSoft) {
  if (fLeptonTagSoftMuons->GetEntries() < 2) return false;
  
  if (fVerbosityLevel > 0) {
    cout << "PhotonTreeWriter::VHLepHasDimuon: Found " 
         << fLeptonTagSoftMuons->GetEntries() << " muons!" << endl;
  }
  
  vector<UInt_t> goodMuons;
  
  // Loop over muons and apply the cut Delta R (mu, pho) > 0.5 .
  for (UInt_t imu=0; imu < fLeptonTagSoftMuons->GetEntries(); ++imu){
    const Muon *mu = fLeptonTagSoftMuons->At(imu);
    if (MathUtils::DeltaR(mu, phHard) < 0.5) continue;
    if (MathUtils::DeltaR(mu, phSoft) < 0.5) continue;
    goodMuons.push_back(imu);
  }
  
  if (goodMuons.size() < 2) return false;

  // Loop over muon pairs of selected muons and apply the cut 70 < mass(mu1, mu2) < 110.
  for (UInt_t iimu1 = 0; iimu1 < goodMuons.size() - 1; ++iimu1) {
    UInt_t imu1 = goodMuons[iimu1];
    const Muon *mu1 = fLeptonTagSoftMuons->At(imu1);
    for (UInt_t iimu2 = iimu1 + 1; iimu2 < goodMuons.size(); ++iimu2) {
      UInt_t imu2 = goodMuons[iimu2];
      const Muon *mu2 = fLeptonTagSoftMuons->At(imu2);
      Double_t mass12 = (mu1->Mom() + mu2->Mom()).M();
      if (mass12 < 70. || 110. < mass12) continue;
      if (fVerbosityLevel > 0) {
        cout << "    Found a tagging dimoun!" << endl << flush;
      }      
      fDiphotonEvent->mu1Pt  = mu1->Pt ();
      fDiphotonEvent->mu1Eta = mu1->Eta();
      fDiphotonEvent->mu1Phi = mu1->Phi();
      fDiphotonEvent->mu2Pt  = mu2->Pt ();
      fDiphotonEvent->mu2Eta = mu2->Eta();
      fDiphotonEvent->mu2Phi = mu2->Phi();
      return true;
    }
  }
  
  return false;
  
} // PhotonTreeWriter::VHLepHasDimuon(..)


//_____________________________________________________________________________
UInt_t PhotonTreeWriter::VHLepNumberOfJets(const Photon *phHard,
                                           const Photon *phSoft,
                                           const Vertex *selvtx,
                                           const Particle *lepton) {

  UInt_t nJets = 0;
  
  // Loop over jets, count those passing selection
  // Use same ID as for the tth tag
  for(UInt_t ijet=0; ijet < fPFJets->GetEntries(); ++ijet){
    const Jet *jet = fPFJets->At(ijet);
    // Apply jet selection, see L116 and L125 of the AN
    if (jet->Pt() < 20. || jet->AbsEta() > 2.4) continue; 
    // Apply the cut Delta R(photon, jet) < 0.5.
    if (MathUtils::DeltaR(jet, phHard) < 0.5) continue;
    if (MathUtils::DeltaR(jet, phSoft) < 0.5) continue;
    // Apply the cut Delta R(photon, lepton) < 0.5.
    if (MathUtils::DeltaR(jet, lepton) < 0.5) continue;
    // Make sure we have a PF jet
    const PFJet *pfjet = dynamic_cast<const PFJet*>(jet);
    if (!pfjet) continue;
    //if (!JetTools::passPFLooseId(pfjet)) continue;
    // Apply the jet ID / pileup removal as given in Table 4
    if (fApplyJetId && !fJetId.passCut(pfjet,selvtx,fPV)) continue;
    // this jet passes, count it in
    ++nJets;
  } // End of loop over jets  
                                                   
  return nJets;
  
} // PhotonTreeWriter::VHLepNumberOfJets(..)


//_____________________________________________________________________________
// Applies the ttH tag given precelected leading and trailing photons
// phHard and phSoft and the corresponding (pre?) selected vertex selvtx. 
// The result is stored as an integer value of the tthTag variable
// entering the diphoton event record.
void PhotonTreeWriter::ApplyTTHTag(const Photon *phHard,
                                   const Photon *phSoft,
                                   const Vertex *selvtx)
{
  // ttH tag = -1 .. not applied
  //            0 .. not tagged
  //            1 .. tagged as a hadronic ttH event
  //            2 .. tagged as a leptonic ttH event
  //            3 .. tagged as both a hadronic and a leptonic ttH event
  
  if (fVerbosityLevel > 90) {
    cout << "+ Entering ApplyTTHTagTag ..." << endl;
    cout << "+ fLeptonTagElectrons->GetEntries(): " << fLeptonTagElectrons->GetEntries() << endl;
    cout << "+ fLeptonTagMuons->GetEntries(): " << fLeptonTagMuons->GetEntries() << endl;
  }
  
  
  //printf("ApplyTTHLepTag: electrons = %i, softelectrons = %i, muons = %i, softmuons = %i\n",fLeptonTagElectrons->GetEntries(), fLeptonTagSoftElectrons->GetEntries(), fLeptonTagMuons->GetEntries(), fLeptonTagSoftMuons->GetEntries());
  
  
  fDiphotonEvent->tthTag = 0;
  
  // Selection taken from the AN2012_480_V6 of 24 April 2013 further
  // refferred to as "the AN"
  // And from "the Hgg AN"
  // http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2013_253_v3.pdf

  const Particle *lepton = TTHSelectLepton(phHard, phSoft, selvtx);

  //printf("lepton = %p\n",(void*)lepton);
  
  // Init jet object counters
  UInt_t nJets = 0;
  UInt_t nBJets = 0;

  // Loop over jets, count those passing selection.
  for(UInt_t ijet=0; ijet < fPFJets->GetEntries(); ++ijet){
    const Jet *jet = fPFJets->At(ijet);
    // Apply jet selection, see L116 and L125 of the AN
    if (jet->Pt() < 25. || jet->AbsEta() > 2.4) continue;
    // Apply Delta R(photon, jet), see email from Francesco Micheli
    // sent 31 Aug 2013
    if (MathUtils::DeltaR(jet, phHard) < 0.5) continue;
    if (MathUtils::DeltaR(jet, phSoft) < 0.5) continue;
    // Apply Delta R(lepton, jet), see emails from Franceso Micheli
    // from 24 and 26 Nov 2013
    if (lepton && MathUtils::DeltaR(jet, lepton) < 0.5) continue;
    // Make sure we have a PF jet
    const PFJet *pfjet = dynamic_cast<const PFJet*>(jet);
    if (!pfjet) continue;
    //if (!JetTools::passPFLooseId(pfjet)) continue;
    // Apply the jet ID as given in Table 4
    if (fApplyJetId && !fJetId.passCut(pfjet,selvtx,fPV)) continue;
    // this jet passes, count it in
    ++nJets;
    // Select b-jets that pass the CSV medium working point, see L128 of the AN
    // and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    if (jet->CombinedSecondaryVertexBJetTagsDisc() < 0.679) continue;
    ++nBJets;
  } // End of loop over jets

  // Check the hadron tag, see Table 7 near L196 of the AN
  if (nJets >= 5 && nBJets >= 1)           fDiphotonEvent->tthTag += 1;

  // Check the lepton tag, see Table 7 near L196 of the AN
  // It has a precedence if both the leptonic and hadronic tags pass.
  // (private e-mail from Francesco Micheli on 15 July 2013).
  if (nJets >= 2 && nBJets >= 1 && lepton) fDiphotonEvent->tthTag += 2;
  
  if (fVerbosityLevel > 90) {
    cout << "+ nJets: " <<  nJets << endl
         << "+ nBJets: " << nBJets << endl
         << "+ lepton: " << lepton << endl
         << "+ tthTag: " << fDiphotonEvent->tthTag << endl;
  }
  
  
} // void PhotonTreeWriter::ApplyTTHTag()


//_____________________________________________________________________________
UInt_t PhotonTreeWriter::NumberOfJets(const Photon *phHard,
                                      const Photon *phSoft,
                                      const Vertex *selvtx,
                                      const double minJetPt,
                                      const double maxAbsEta) {

  UInt_t nJets = 0;

  // Loop over jets, count those passing selection
  // Use same ID as for the tth tag
  for(UInt_t ijet=0; ijet < fPFJets->GetEntries(); ++ijet){
    const Jet *jet = fPFJets->At(ijet);
    // Apply jet selection, see L116 and L125 of the AN
    if (jet->Pt() < minJetPt || jet->AbsEta() > maxAbsEta) continue;
    // Apply the cut Delta R(photon, jet) < 0.5.
    if (MathUtils::DeltaR(jet, phHard) < 0.5) continue;
    if (MathUtils::DeltaR(jet, phSoft) < 0.5) continue;
    // Make sure we have a PF jet
    const PFJet *pfjet = dynamic_cast<const PFJet*>(jet);
    if (!pfjet) continue;
    //if (!JetTools::passPFLooseId(pfjet)) continue;
    // Apply the jet ID / pileup removal as given in Table 4
    if (fApplyJetId && !fJetId.passCut(pfjet,selvtx,fPV)) continue;
    // this jet passes, count it in
    ++nJets;
  } // End of loop over jets

  return nJets;

} // NumberOfJets


//_____________________________________________________________________________
// PhotonTreeWriter::PFJetVector *
// PhotonTreeWriter::GetSelectedPFJets(const DeltaRVetoVector &drVetos,
//                                     const Vertex &vertex,
//                                     const double minPt,
//                                     const double maxAbsEta) {
// 
//   PFJetVector *pfjets = new PFJetVector();
// 
//   // Loop over jets, count those passing selection
//   // Use same ID as for the tth tag
//   for(UInt_t ijet=0; ijet < fPFJets->GetEntries(); ++ijet){
//     const Jet *jet = fPFJets->At(ijet);
//     // Apply jet selection, see L116 and L125 of the AN
//     if (jet->Pt() < minPt || jet->AbsEta() > maxAbsEta) continue;
//     // Apply the vetos Delta R(jet, particle) > minDeltaR.
//     DeltaRVetoVector::const_iterator drVeto = drVetos.begin();
//     for (; drVeto < drVetos.end(); ++drVeto) {
//       const Particle *particle = drVeto->first;
//       std::cout << "run: " << fDiphotonEvent->run
//                 << ", lumi: " << fDiphotonEvent->lumi
//                 << ", events: " << fDiphotonEvent->evt
//                 /// The line below causes a segmentation fault!
//                 << ", veto particle pt: " << particle->Mom().Pt()
//                 << ", eta: " << particle->Mom().Eta()
//                 << ", phi: " << particle->Mom().Phi()
//                 << "\n";
//       double minDeltaR = drVeto->second;
//       if (MathUtils::DeltaR(particle, jet) < minDeltaR) break;
//     } /// Loop over Delta R vetos.
//     if (drVeto < drVetos.end()) continue; /// failed a Delta R veto
//     // Make sure we have a PF jet
//     const PFJet *pfjet = dynamic_cast<const PFJet*>(jet);
//     if (!pfjet) continue;
//     if (!JetTools::passPFLooseId(pfjet)) continue;
//     // Apply the jet ID / pileup removal as given in Table 4
//     Double_t betaStar = JetTools::betaStarClassic(pfjet, &vertex, fPV);
//     if (betaStar > 0.2 * log(fPV->GetEntries() - 0.64)) continue;
//     if (JetTools::dR2Mean(pfjet, -1) > 0.065) continue;
//     // this jet passes, count it in
//     pfjets->push_back(pfjet);
//   } // End of loop over jets
// 
//   return pfjets;
// 
// } // GetSelectedPFJets


//_____________________________________________________________________________
UInt_t PhotonTreeWriter::NumberOfBJets(const Photon *phHard,
                                       const Photon *phSoft,
                                       const Vertex *selvtx,
                                       const double minJetPt,
                                       const double maxAbsEta) {

  UInt_t nBJets = 0;

  // Loop over jets, count those passing selection
  // Use same ID as for the tth tag
  for(UInt_t ijet=0; ijet < fPFJets->GetEntries(); ++ijet){
    const Jet *jet = fPFJets->At(ijet);
    // Apply jet selection, see L116 and L125 of the AN
    if (jet->Pt() < minJetPt || jet->AbsEta() > maxAbsEta) continue;
    // Apply the cut Delta R(photon, jet) < 0.5.
    if (MathUtils::DeltaR(jet, phHard) < 0.5) continue;
    if (MathUtils::DeltaR(jet, phSoft) < 0.5) continue;
    // Make sure we have a PF jet
    const PFJet *pfjet = dynamic_cast<const PFJet*>(jet);
    if (!pfjet) continue;
    //if (!JetTools::passPFLooseId(pfjet)) continue;
    // Apply the jet ID / pileup removal as given in Table 4
    if (fApplyJetId && !fJetId.passCut(pfjet,selvtx,fPV)) continue;
    // Select b-jets that pass the CSV medium working point, see L128 of the AN
    // and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    if (jet->CombinedSecondaryVertexBJetTagsDisc() < 0.679) continue;
    // this jet passes, count it in
    ++nBJets;
    if (fDoSynching && fDiphotonEvent->bjetcsv <= -99) {
      fDiphotonEvent->bjetcsv = jet->CombinedSecondaryVertexBJetTagsDisc();
    }
  } // End of loop over jets

  return nBJets;
} // NumberOfBJets


//_____________________________________________________________________________
const Particle *
PhotonTreeWriter::TTHSelectLepton(const Photon *phHard,
                                  const Photon *phSoft,
                                  const Vertex *selvtx)
{
  if (fVerbosityLevel > 90) {
    cout << "++ Entering TTHSelectLepton ..." << endl;
  }   
  const Electron * electron = TTHSelectElectron(phHard, phSoft, selvtx);
  const Muon *     muon     = TTHSelectMuon    (phHard, phSoft, selvtx);
  const Particle * lepton = 0;
  /// Select the lepton with the higher Pt (private discussion
  /// with Francesco Micheli on 26 Nov 2013)
  if (electron && muon) {
    if (electron->Pt() > muon->Pt()) {lepton = electron;}
    else                             {lepton = muon    ;}
  } else if (electron && !muon) {
    lepton = electron;
  } else if (!electron && muon) {
    lepton = muon    ;
  }
  if (fVerbosityLevel > 90) {
    cout << "++ electron: " << electron << endl
         << "++ muon: " << muon << endl 
         << "++ lepton: " << lepton << endl
         << "++ Exiting TTHSelectLepton ..." << endl;
  }
  return lepton;
} // TTHSelectLepton


//_____________________________________________________________________________
// Returns the highest-MVA-ID electron passing all cuts, 0 if none found.
const Electron *
PhotonTreeWriter::TTHSelectElectron(const Photon *phHard,
                                    const Photon *phSoft,
                                    const Vertex *selvtx)
{
  if (fVerbosityLevel > 90) {
    cout << "+++ Entering TTHSelectElectron ..." << endl
         << "+++ PhotonTools::ElectronVetoCiC(phHard, fLeptonTagElectrons): " 
         <<      PhotonTools::ElectronVetoCiC(phHard, fLeptonTagElectrons) << endl
         << "+++ PhotonTools::ElectronVetoCiC(phSoft, fLeptonTagElectrons): " 
         <<      PhotonTools::ElectronVetoCiC(phSoft, fLeptonTagElectrons) << endl
         << "+++ PhotonTools::ElectronVetoCiC(phHard, fElectrons): " 
         <<      PhotonTools::ElectronVetoCiC(phHard, fElectrons) << endl
         << "+++ PhotonTools::ElectronVetoCiC(phSoft, fElectrons): " 
         <<      PhotonTools::ElectronVetoCiC(phSoft, fElectrons) << endl;
  }   
  
  if (PhotonTools::ElectronVetoCiC(phHard, fElectrons) >= 1.          &&
      PhotonTools::ElectronVetoCiC(phSoft, fElectrons) >= 1.){
    // Loop over electrons, apply all cuts, find the one wiht hightes ID MVA
    for (UInt_t iele=0; iele < fLeptonTagElectrons->GetEntries(); ++iele) {
      const Electron *ele = fLeptonTagElectrons->At(iele);

      if (fVerbosityLevel > 90) {
        cout << "+++ iele: " << iele << endl
             << "+++ ele: " << ele << endl
             << "+++ ele->Pt(): " << ele->Pt() << endl
             << "+++ ele->AbsEta(): " << ele->AbsEta() << endl
             << "+++ MathUtils::DeltaR(ele, phHard): " 
             <<      MathUtils::DeltaR(ele, phHard) << endl
             << "+++ MathUtils::DeltaR(ele, phSoft): " 
             <<      MathUtils::DeltaR(ele, phSoft) << endl
             << "+++ MassOfPairIsWithinWindowAroundMZ(ele, phHard, 10): " 
             <<      MassOfPairIsWithinWindowAroundMZ(ele, phHard, 10) << endl
             << "+++ MassOfPairIsWithinWindowAroundMZ(ele, phSoft, 10): " 
             <<      MassOfPairIsWithinWindowAroundMZ(ele, phSoft, 10) << endl;
             //<< "+++ GetElectronIdMva(ele): " << GetElectronIdMva(ele) << endl;
      }   
      
      // Apply kinematic cuts, see L133 and L134 of the AN
      //(JOSH: In fact these are already applied in the ElectronIDMod as we have configured it)
      //if (ele->Pt() < 20. || ele->AbsEta() > 2.5) continue;
      //keep Pt cut here in any case we reorganize logic later
      if (ele->Pt()<20.) continue;
      
      // Require separation between this electron and both photons,
      // see the slide 7, bullet 2
      if (MathUtils::DeltaR(ele, phHard) < 1.) continue;
      if (MathUtils::DeltaR(ele, phSoft) < 1.) continue;
      // Require electron-photon mass outside of a 20 GeV window around MZ
      //JOSH: Shouldn't be applied for ttH category
//       if (MassOfPairIsWithinWindowAroundMZ(ele, phHard, 10)) continue;
//       if (MassOfPairIsWithinWindowAroundMZ(ele, phSoft, 10)) continue;
      // Electron ID MVA arbitration (private discussion with
      // Francesco Micheli on 26 Nov 2013)
//       double idMva = GetElectronIdMva(ele);
//       if (idMva > maxIdMva) {
//         maxIdMva = idMva;
//         selectedElectron = ele;
//       }
      //JOSH: electrons here are already sorted by mva output, so can just take the first one and break
      //(in fact the call to GetElectronIdMva was segfaulting for a reason which I didn't care to debug)
      return ele;
    } // Loop over electrons
  }

  return 0;
  
} // TTHSelectElectron


//_____________________________________________________________________________
// Returns the highest-pt muon passing all the cuts, 0 if none found.
const Muon *
PhotonTreeWriter::TTHSelectMuon(const Photon *phHard,
                                const Photon *phSoft,
                                const Vertex *selvtx)
{
  if (fVerbosityLevel > 90) {
    cout << "+++ Entering TTHSelectMuon ..." << endl;
  }   
  const Muon *selectedMuon = 0;
  // Loop over muons
  for (UInt_t imu=0; imu < fLeptonTagMuons->GetEntries(); ++imu) {
    const Muon *mu = fLeptonTagMuons->At(imu);
    // Apply kinematic cuts, see L132 and L134 of the AN
    if (mu->Pt() < 20. || mu->AbsEta() > 2.4) continue;
    // Same as for electrons, require separation between this muon and both
    // photons.
    // Also confirmed by Francesco Micheli in an e-mail from 15 July 2013
    if (MathUtils::DeltaR(mu, phHard) < 0.5) continue;
    if (MathUtils::DeltaR(mu, phSoft) < 0.5) continue;
    // The first muon found is the highest-pt one because muons are pt-ordered.
    selectedMuon = mu;
    break;
  }
  if (fVerbosityLevel > 90) {
    cout << "+++ fLeptonTagMuons->GetEntries(): " 
         << fLeptonTagMuons->GetEntries() << endl
         << "+++ Looping over muons ..." << endl;
    for (UInt_t imu=0; imu < fLeptonTagMuons->GetEntries(); ++imu) {
      const Muon *mu = fLeptonTagMuons->At(imu);
      cout << "++++ imu: " << imu << endl
           << "++++ mu: " << mu << endl
           << "++++ mu pt: " << mu->Pt() << endl
           << "++++ mu |eta|: " << mu->AbsEta() << endl
           << "++++ MathUtils::DeltaR(mu, phHard): " 
           << MathUtils::DeltaR(mu, phHard) << endl
           << "++++ MathUtils::DeltaR(mu, phSoft): " 
           << MathUtils::DeltaR(mu, phSoft) << endl;
    }    
    cout << "+++ Finished loop over muons ..." << endl
         << "+++ selectedMuon: " << selectedMuon << endl
         << "+++ Exiting TTHSelectMuon ..." << endl << flush;
  }   
  return selectedMuon;
} // TTHSelectMuon


//_____________________________________________________________________________
template<typename Element, typename Collection>
UInt_t PhotonTreeWriter::IndexOfElementInCollection(
          const Element * element,
          const Collection * collection
          )
{
  UInt_t index = 0;
  for (; index < collection->GetEntries() && index < (UInt_t) -1; ++index) {
    if (element == collection->At(index)) {
      break;
    }
  }
  return index;
}

//_____________________________________________________________________________
UInt_t PhotonTreeWriter::IndexOfNearestSuperClusterInCollection(
          const SuperCluster    *element,
          const SuperClusterCol *collection
          )
{
  double minMass = 999.;
  UInt_t minIndex = 0;
  if (collection->GetEntries() > 0) {
    minMass = SuperClusterPairMass(element, collection->At(0));
  }
  for (UInt_t index = 1;
       index < collection->GetEntries() && index < (UInt_t) -1; ++index) {
    double mass = SuperClusterPairMass(element, collection->At(index));
    if (mass < minMass) {
      minMass = mass;
      minIndex = index;
    }
  }
  return minIndex;
}


//_____________________________________________________________________________
double PhotonTreeWriter::SuperClusterPairMass(const SuperCluster *sc1,
                                              const SuperCluster *sc2)
{
  FourVectorM p1 = SuperClusterFourVectorM(sc1);
  FourVectorM p2 = SuperClusterFourVectorM(sc2);
  return (p1 + p2).M();
}


//_____________________________________________________________________________
FourVectorM PhotonTreeWriter::SuperClusterFourVectorM(const SuperCluster *sc)
{
  double e   = sc->Energy();
  double eta = sc->Eta();
  double phi = sc->Phi();
  double pt  = e / TMath::CosH(eta);
  return FourVectorM(pt, eta, phi, 0.);
}


//_____________________________________________________________________________
void PhotonTreeWriter::Terminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}


//_____________________________________________________________________________
void PhotonTreeWriter::LogEventInfo()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, 80, "%c", timeinfo);

  const EventHeader *event = GetEventHeader();
  cout << buffer << ": Processing " << fProcessedEvents << ". record"
       << ", run " <<  event->RunNum()
       << ", lumi " << event->LumiSec()
       << ", event " << event->EvtNum() << endl;  
}


//______________________________________________________________________________
double
PhotonTreeWriter::GetElectronIdMva(const Electron *electron)
{
  Int_t closestVtx = 0;
  Double_t distVtx = 999.0;
  /// Find the closest vertex in dz.
  for(UInt_t nv=0; nv<fPV->GetEntries(); nv++){
    double dz = TMath::Abs(electron->GsfTrk()->DzCorrected(*fPV->At(nv)));
    if(dz < distVtx) {
      distVtx    = dz;
      closestVtx = nv;
    }
  }
  return fElectronIDMVA->MVAValue(electron, fPV->At(closestVtx));
} // GetElectronIdMva
