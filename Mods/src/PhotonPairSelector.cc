#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "TDataMember.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>

using namespace mithep;

ClassImp(mithep::PhotonPairSelector)
ClassImp(mithep::PhotonPairSelectorPhoton)
ClassImp(mithep::PhotonPairSelectorDiphotonEvent)


//--------------------------------------------------------------------------------------------------
PhotonPairSelector::PhotonPairSelector(const char *name, const char *title) : 
  // Base Module...
  BaseMod            (name,title),

  // define all the Branches to load
  fPhotonBranchName  (Names::gkPhotonBrn),
  fElectronName      (Names::gkElectronBrn),
  fConversionName    (Names::gkMvfConversionBrn),  
  fTrackBranchName   (Names::gkTrackBrn),
  fPileUpDenName     (Names::gkPileupEnergyDensityBrn),
  fPVName            (Names::gkPVBeamSpotBrn),
  fBeamspotName      (Names::gkBeamSpotBrn),
  fPFCandName        (Names::gkPFCandidatesBrn),
  // MC specific stuff...
  fMCParticleName    (Names::gkMCPartBrn),
  fPileUpName        (Names::gkPileupInfoBrn),

  fGoodPhotonsName   (ModNames::gkGoodPhotonsName),

  // ----------------------------------------
  // Selection Types
  fPhotonSelType     ("NoSelection"),
  fVertexSelType     ("StdSelection"),
  fPhSelType         (kNoPhSelection),
  fVtxSelType        (kStdVtxSelection),
  
  // ----------------------------------------
  fPhotonPtMin       (20.0),
  fPhotonEtaMax      (2.5),

  fLeadingPtMin      (40.0),
  fTrailingPtMin     (30.0),

  fIsData            (false),
  fPhotonsFromBranch (true),  
  fPVFromBranch      (true),
  
  // ----------------------------------------
  // collections....
  fPhotons           (0),
  fElectrons         (0),
  fConversions       (0),
  fTracks            (0),
  fPileUpDen         (0),
  fPV                (0),
  fBeamspot          (0),
  fPFCands           (0),
  fMCParticles       (0),
  fPileUp            (0),

  // ---------------------------------------
  fDataEnCorr_EB_hR9 (0.),
  fDataEnCorr_EB_lR9 (0.),
  fDataEnCorr_EE_hR9 (0.),
  fDataEnCorr_EE_lR9 (0.),

  fRunStart          (0),
  fRunEnd            (0),

  fMCSmear_EB_hR9    (0.),
  fMCSmear_EB_lR9    (0.),
  fMCSmear_EE_hR9    (0.),
  fMCSmear_EE_lR9    (0.),

  // ---------------------------------------
  rng                (new TRandom3()),  

  // ---------------------------------------
  fDoDataEneCorr     (true),
  fDoMCSmear         (true),
  fDoVtxSelection    (true),
  fApplyEleVeto      (true),

  fTupleName         ("hCiCtuple")
  
{
  // Constructor.
}

PhotonPairSelector::~PhotonPairSelector(){
  if(rng) delete rng;
}

//--------------------------------------------------------------------------------------------------
void PhotonPairSelector::Process()
{
  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fConversionName,     fConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);    
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);



  // -----------------------------------------------------------
  // OUtput Photon Collection. It will ALWAYS conatrin either 0 or 2 Photons
  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(kTRUE);
  

  // ------------------------------------------------------------  
  // load event based information
  Int_t _numPU = -1.;        // some sensible default values....
  Int_t _numPUminus = -1.;        // some sensible default values....
  Int_t _numPUplus = -1.;        // some sensible default values....

  Float_t _tRho  = -99.;
  if( fPileUpDen->GetEntries() > 0 )
    _tRho  = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();


  
  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));

  // ------------------------------------------------------------  
  // Get Event header for Run info etc.
  const EventHeader* evtHead = this->GetEventHeader();
  unsigned int evtNum = evtHead->EvtNum();
  Float_t _evtNum1   = (Float_t) ( (int) (evtNum/10000.) );
  Float_t _evtNum2   = (Float_t) ( (int) (evtNum % 10000)  );
  UInt_t   runNumber = evtHead->RunNum();
  Float_t _runNum    = (Float_t) runNumber;
  Float_t _lumiSec   = (Float_t) evtHead->LumiSec();
  
  // ------------------------------------------------------------
  // here we'll store the preselected Photons (and which CiCCategory they are...)
  PhotonOArr* preselPh  = new PhotonOArr;
  std::vector<PhotonTools::CiCBaseLineCats> preselCat;
  preselCat.resize(0);
  
  // 1. we do the pre-selection; but keep the non-passing photons in a secont container...
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);

    if(ph->SCluster()->AbsEta()>= fPhotonEtaMax || (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566)) continue;
    if(ph->Et()                <  fPhotonPtMin)     continue;
    if(ph->HadOverEm()         >  0.15)     continue;
    if(ph->IsEB()) {
      if(ph->CoviEtaiEta() > 0.013) continue;      
    } else {
      if(ph->CoviEtaiEta() > 0.03) continue;
    }    
    preselPh->Add(ph);
  }

  // Sorry... need the second loop here in order to sort & assign the right Categories..
  preselPh->Sort();
  for(unsigned int iPh = 0; iPh <preselPh->GetEntries(); ++iPh) {
    const Photon* ph = preselPh->At(iPh);
    preselCat.push_back(PhotonTools::CiCBaseLineCat(ph));
  }
  
  // ------------------------------------------------------------
  // compute how many pairs there are ...
  unsigned int numPairs = 0;
  if( preselPh->GetEntries() > 0) numPairs = (preselPh->GetEntries()-1)*preselPh->GetEntries()/2;  
  // ... and create all possible pairs of pre-selected photons
  std::vector<unsigned int> idx1st;
  std::vector<unsigned int> idx2nd;
  std::vector<PhotonTools::CiCBaseLineCats> cat1st;
  std::vector<PhotonTools::CiCBaseLineCats> cat2nd;
  // ... this will be used to store whether a givne pair passes the cuts
  std::vector<bool> pairPasses;
  
  if(numPairs > 0) {
    for(unsigned int i1st = 0; i1st <preselPh->GetEntries() - 1; ++i1st) {
      for(unsigned int i2nd = i1st + 1; i2nd <preselPh->GetEntries(); ++i2nd) {
	idx1st.push_back(i1st);
	idx2nd.push_back(i2nd);
	pairPasses.push_back(true);
      }
    }
  }

  // ------------------------------------------------------------  
  // array to store the index of 'chosen Vtx' for each pair
  const Vertex** theVtx        = new const Vertex*[numPairs];    // holds the 'chosen' Vtx for each Pair
  Photon**       fixPh1st      = new       Photon*[numPairs];    // holds the 1st Photon for each Pair       
  Photon**       fixPh2nd      = new       Photon*[numPairs];    // holds the 2nd photon for each Pair
  

  // store pair-indices for pairs passing the selection
  std::vector<unsigned int> passPairs;
  passPairs.resize(0);
  
  // ------------------------------------------------------------  
  // Loop over all Pairs and to the 'incredible machine' running....
  for(unsigned int iPair = 0; iPair < numPairs; ++iPair) {    

    // first we need a hard copy of the incoming photons
    fixPh1st[iPair] = new Photon(*preselPh->At(idx1st[iPair]));
    fixPh2nd[iPair] = new Photon(*preselPh->At(idx2nd[iPair]));
    // we also store the category, so we don't have to ask all the time...
    cat1st.push_back(preselCat[idx1st[iPair]]);
    cat2nd.push_back(preselCat[idx2nd[iPair]]);
    
    // now we dicide if we either scale (Data) or Smear (MC) the Photons
    if (fIsData) {
      if(fDoDataEneCorr) {
	// statring with scale = 1.
	double scaleFac1 = 1.;
	double scaleFac2 = 1.;
	// checking the run Rangees ...
	Int_t runRange = FindRunRangeIdx(runNumber);
	if(runRange > -1) { 
	  scaleFac1 += GetDataEnCorr(runRange, cat1st[iPair]);
	  scaleFac2 += GetDataEnCorr(runRange, cat2nd[iPair]);
	}      
	PhotonTools::ScalePhoton(fixPh1st[iPair], scaleFac1);
	PhotonTools::ScalePhoton(fixPh2nd[iPair], scaleFac2);
      }
    } else {
      if(fDoMCSmear) {      
	// get the seed to do deterministic smearing...
	UInt_t seedBase = (UInt_t) evtNum + (UInt_t) _runNum + (UInt_t) _lumiSec;
	UInt_t seed1    = seedBase + (UInt_t) fixPh1st[iPair]->E() + (UInt_t) (TMath::Abs(10.*fixPh1st[iPair]->SCluster()->Eta()));
	UInt_t seed2    = seedBase + (UInt_t) fixPh2nd[iPair]->E() + (UInt_t) (TMath::Abs(10.*fixPh2nd[iPair]->SCluster()->Eta()));
	// get the smearing for MC photons..

	double width1 = GetMCSmearFac(cat1st[iPair]);
	double width2 = GetMCSmearFac(cat2nd[iPair]);
	
	PhotonTools::SmearPhoton(fixPh1st[iPair], rng, width1, seed1);
	PhotonTools::SmearPhoton(fixPh2nd[iPair], rng, width2, seed2);

      }
    }

    // store the vertex for this pair
    switch( fVtxSelType ){
    case kStdVtxSelection:
      theVtx[iPair] = fPV->At(0);
      break;
    case kCiCVtxSelection:
      theVtx[iPair] = VertexTools::findVtxBasicRanking(fixPh1st[iPair],fixPh2nd[iPair], bsp, fPV, fConversions);
      break;
    case kMITVtxSelection:
      // need PFCandidate Collection
      theVtx[iPair] = VertexTools::BestVtx(fPFCands, fPV, bsp, mithep::FourVector((fixPh1st[iPair]->Mom()+fixPh2nd[iPair]->Mom()))); 
      break;
    default:
      theVtx[iPair] = fPV->At(0);

    }

    // fix the kinematics for both events
    FourVectorM newMom1st = fixPh1st[iPair]->MomVtx(theVtx[iPair]->Position());
    FourVectorM newMom2nd = fixPh2nd[iPair]->MomVtx(theVtx[iPair]->Position());
    fixPh1st[iPair]->SetMom(newMom1st.X(), newMom1st.Y(), newMom1st.Z(), newMom1st.E());
    fixPh2nd[iPair]->SetMom(newMom2nd.X(), newMom2nd.Y(), newMom2nd.Z(), newMom2nd.E());

    // check if both photons pass the CiC selection
    // FIX-ME: Add other possibilities....
    bool pass1 = false;
    bool pass2 = false;
    switch( fVtxSelType ){
    case kNoPhSelection:
      pass1 = ( fixPh1st[iPair]->Pt() > fLeadingPtMin  );
      pass2 = ( fixPh2nd[iPair]->Pt() > fTrailingPtMin );
      break;
    case kCiCPhSelection:


      pass1 = PhotonTools::PassCiCSelection(fixPh1st[iPair], theVtx[iPair], fTracks, fElectrons, fPV, _tRho, fLeadingPtMin, fApplyEleVeto);
      if(pass1) pass2 = PhotonTools::PassCiCSelection(fixPh2nd[iPair], theVtx[iPair], fTracks, fElectrons, fPV, _tRho, fTrailingPtMin, fApplyEleVeto);

      break;
    case kMITPhSelection:
      // FIX-ME: This is a place-holder.. MIT guys: Please work hard... ;)
      pass1 = ( fixPh1st[iPair]->Pt() > fLeadingPtMin  );
      pass2 = ( fixPh2nd[iPair]->Pt() > fTrailingPtMin );
      break;
    default:
      pass1 = true;
      pass2 = true;
    }
    // finally, if both Photons pass the selections, add the pair to the 'passing Pairs)
    if( pass1 && pass2 ) passPairs.push_back(iPair);
  }
  
  
  // ---------------------------------------------------------------
  // ... we're almost done, stau focused...
  // loop over all passing pairs and find the one with the highest sum Et
  const Vertex* _theVtx  = NULL;
  Photon*        phHard  = NULL;
  Photon*        phSoft  = NULL;

  PhotonTools::CiCBaseLineCats catPh1 = PhotonTools::kCiCNoCat;
  PhotonTools::CiCBaseLineCats catPh2 = PhotonTools::kCiCNoCat;
  
  double maxSumEt = 0.;
  for(unsigned int iPair=0; iPair<passPairs.size(); ++iPair){
    double sumEt = fixPh1st[passPairs[iPair]]->Et();
    sumEt += fixPh2nd[passPairs[iPair]]->Et();
    if( sumEt > maxSumEt ) {
      maxSumEt = sumEt;
      phHard = fixPh1st[passPairs[iPair]];
      phSoft = fixPh2nd[passPairs[iPair]];
      catPh1 = cat1st[passPairs[iPair]];
      catPh2 = cat2nd[passPairs[iPair]];
      _theVtx = theVtx[iPair];
    }
  }
  
  // ---------------------------------------------------------------
  // we have the Photons (*PARTY*)... compute some useful qunatities

  Float_t _theVtxZ = -999.;

  if(phHard && phSoft) {
    GoodPhotons->AddOwned(phHard);
    GoodPhotons->AddOwned(phSoft);
  }

  if(_theVtx) _theVtxZ=_theVtx->Position().Z();
  
  // sort according to pt
  GoodPhotons->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);

  // delete auxiliary photon collection...
  delete preselPh;
  delete[] theVtx;
  
  // Fill the useful information (mass and di-photon pt)
  bool doFill   = ( phHard && phSoft );
  if (!doFill) return;

  if( !fIsData ) {
    LoadBranch(fMCParticleName);
    LoadBranch(fPileUpName);
  }  
  
  if( !fIsData ) {
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing()==0) _numPU = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() == -1) _numPUminus = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() ==  1) _numPUplus  = puinfo->GetPU_NumInteractions();
    }
  }

  Float_t _mass = ( doFill ? (phHard->Mom()+phSoft->Mom()).M()  : -100.);
  Float_t _ptgg = ( doFill ? (phHard->Mom()+phSoft->Mom()).Pt() : -100.);
  
  // in case of a MC event, try to find Higgs and Higgs decay Z poisition
  Float_t _pth    = -100.;
  Float_t _decayZ = -100.;
  if( !fIsData ) FindHiggsPtAndZ(_pth, _decayZ);
  
  Float_t evtCat = -1.;
  if( doFill ) {
    evtCat = GetEventCat(catPh1, catPh2);
    if(_ptgg < 40.) evtCat += 4.;
  }
  
  const MCParticle *phgen1 = PhotonTools::MatchMC(phHard,fMCParticles);
  const MCParticle *phgen2 = PhotonTools::MatchMC(phSoft,fMCParticles);

  Float_t _gencostheta = -99.;
  if (phgen1 && phgen2) {
    _gencostheta = ThreeVector(phgen1->Mom()).Unit().Dot(ThreeVector(phgen2->Mom()).Unit());
  }
  
   if(_mass > 0.) {
    fDiphotonEvent->rho = _tRho;
    fDiphotonEvent->genHiggspt = _pth;
    fDiphotonEvent->genHiggsZ = _decayZ;
    fDiphotonEvent->gencostheta = _gencostheta;
    fDiphotonEvent->vtxZ = _theVtxZ;
    fDiphotonEvent->nVtx = fPV->GetEntries();
    fDiphotonEvent->numPU = _numPU;
    fDiphotonEvent->numPUminus = _numPUminus;
    fDiphotonEvent->numPUplus = _numPUplus;
    fDiphotonEvent->mass = _mass;
    fDiphotonEvent->ptgg = _ptgg;
    fDiphotonEvent->costheta =  ThreeVector(phHard->Mom()).Unit().Dot(ThreeVector(phSoft->Mom()).Unit());
    fDiphotonEvent->evt = GetEventHeader()->EvtNum();
    fDiphotonEvent->run = GetEventHeader()->RunNum();
    fDiphotonEvent->lumi = GetEventHeader()->LumiSec();
    fDiphotonEvent->evtcat = evtCat;

    fDiphotonEvent->photons[0].SetVars(phHard,phgen1);
    fDiphotonEvent->photons[1].SetVars(phSoft,phgen2);
    
    Float_t ph1ecor    = fDiphotonEvent->photons[0].Ecor();
    Float_t ph1ecorerr = fDiphotonEvent->photons[0].Ecorerr();
    Float_t ph2ecor    = fDiphotonEvent->photons[1].Ecor();
    Float_t ph2ecorerr = fDiphotonEvent->photons[1].Ecorerr();
    
    fDiphotonEvent->masscor = TMath::Sqrt(2.0*ph1ecor*ph2ecor*(1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscorerr = 0.5*fDiphotonEvent->masscor*TMath::Sqrt(ph1ecorerr*ph1ecorerr/ph1ecor/ph1ecor + ph2ecorerr*ph2ecorerr/ph2ecor/ph2ecor);
    
    //printf("r9 = %5f, photon sigieie = %5f, seed sigieie = %5f\n",phHard->R9(),phHard->CoviEtaiEta(),sqrt(phHard->SCluster()->Seed()->CoviEtaiEta()));
    
    hCiCTuple->Fill();  
    
    fSinglePhoton->SetVars(phHard,phgen1);
    hCiCTupleSingle->Fill();  
    fSinglePhoton->SetVars(phSoft,phgen2);
    hCiCTupleSingle->Fill();  
    
   }    
    
  return;

}

//--------------------------------------------------------------------------------------------------
void PhotonPairSelector::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  ReqEventObject(fPhotonBranchName,   fPhotons,    fPhotonsFromBranch);
  ReqEventObject(fTrackBranchName,    fTracks,     true);
  ReqEventObject(fElectronName,       fElectrons,  true);  
  ReqEventObject(fPileUpDenName,      fPileUpDen,  true);
  ReqEventObject(fPVName,             fPV,         fPVFromBranch);
  ReqEventObject(fConversionName,     fConversions,true);
  ReqEventObject(fBeamspotName,       fBeamspot,   true);
  ReqEventObject(fPFCandName,         fPFCands,    true);
  
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }
  
  if      (fPhotonSelType.CompareTo("CiCSelection") == 0) 
    fPhSelType =       kCiCPhSelection;
  else if (fPhotonSelType.CompareTo("MITSelection") == 0) 
    fPhSelType =       kMITPhSelection;
  else 
    fPhSelType =       kNoPhSelection;

  if      (fVertexSelType.CompareTo("CiCSelection") == 0) 
    fVtxSelType =       kCiCVtxSelection;
  else if (fVertexSelType.CompareTo("MITSelection") == 0) 
    fVtxSelType =       kMITVtxSelection;
  else 
    fVtxSelType =       kStdVtxSelection;  

  
  //initialize photon energy corrections
  //PhotonFix::initialise("4_2",std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFix.dat")).Data()));  

  fDiphotonEvent = new PhotonPairSelectorDiphotonEvent;
  fSinglePhoton = new PhotonPairSelectorPhoton;

  
  hCiCTuple = new TTree(fTupleName.Data(),fTupleName.Data());
  TString singlename = fTupleName + TString("Single");
  hCiCTupleSingle = new TTree(singlename,singlename);
  
  //make flattish tree from classes so we don't have to rely on dictionaries for reading later
  TClass *eclass = TClass::GetClass("mithep::PhotonPairSelectorDiphotonEvent");
  TClass *pclass = TClass::GetClass("mithep::PhotonPairSelectorPhoton");
  TList *elist = eclass->GetListOfDataMembers();
  TList *plist = pclass->GetListOfDataMembers();
    
  for (int i=0; i<elist->GetEntries(); ++i) {
    const TDataMember *tdm = static_cast<const TDataMember*>(elist->At(i));
    if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
    TString typestring;
    if (TString(tdm->GetTypeName())=="Char_t") typestring = "B";
    else if (TString(tdm->GetTypeName())=="UChar_t") typestring = "b";
    else if (TString(tdm->GetTypeName())=="Short_t") typestring = "S";
    else if (TString(tdm->GetTypeName())=="UShort_t") typestring = "s";
    else if (TString(tdm->GetTypeName())=="Int_t") typestring = "I";
    else if (TString(tdm->GetTypeName())=="UInt_t") typestring = "i";
    else if (TString(tdm->GetTypeName())=="Float_t") typestring = "F";
    else if (TString(tdm->GetTypeName())=="Double_t") typestring = "D";
    else if (TString(tdm->GetTypeName())=="Long64_t") typestring = "L";
    else if (TString(tdm->GetTypeName())=="ULong64_t") typestring = "l";
    else if (TString(tdm->GetTypeName())=="Bool_t") typestring = "O";
    else continue;
    //printf("%s %s: %i\n",tdm->GetTypeName(),tdm->GetName(),int(tdm->GetOffset()));
    Char_t *addr = (Char_t*)fDiphotonEvent;
    assert(sizeof(Char_t)==1);
    hCiCTuple->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
    hCiCTupleSingle->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
  }

  for (int iph=0; iph<2; ++iph) {
    for (int i=0; i<plist->GetEntries(); ++i) {
      const TDataMember *tdm = static_cast<const TDataMember*>(plist->At(i));
      if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
      TString typestring;
      if (TString(tdm->GetTypeName())=="Char_t") typestring = "B";
      else if (TString(tdm->GetTypeName())=="UChar_t") typestring = "b";
      else if (TString(tdm->GetTypeName())=="Short_t") typestring = "S";
      else if (TString(tdm->GetTypeName())=="UShort_t") typestring = "s";
      else if (TString(tdm->GetTypeName())=="Int_t") typestring = "I";
      else if (TString(tdm->GetTypeName())=="UInt_t") typestring = "i";
      else if (TString(tdm->GetTypeName())=="Float_t") typestring = "F";
      else if (TString(tdm->GetTypeName())=="Double_t") typestring = "D";
      else if (TString(tdm->GetTypeName())=="Long64_t") typestring = "L";
      else if (TString(tdm->GetTypeName())=="ULong64_t") typestring = "l";
      else if (TString(tdm->GetTypeName())=="Bool_t") typestring = "O";
      else continue;
      //printf("%s\n",tdm->GetTypeName());
      TString varname = TString::Format("ph%d.%s",iph+1,tdm->GetName());
      
      Char_t *addr = (Char_t*)&fDiphotonEvent->photons[iph];
      assert(sizeof(Char_t)==1);
      hCiCTuple->Branch(varname,addr+tdm->GetOffset(),TString::Format("%s/%s",varname.Data(),typestring.Data()));
      
      if (iph==0) {
        TString singlename = TString::Format("ph.%s",tdm->GetName());
        Char_t *addrsingle = (Char_t*)fSinglePhoton;
        hCiCTupleSingle->Branch(singlename,addrsingle+tdm->GetOffset(),TString::Format("%s/%s",singlename.Data(),typestring.Data()));
      }
    }
  }


  AddOutput(hCiCTuple);
  AddOutput(hCiCTupleSingle);


}

// ----------------------------------------------------------------------------------------
// some helpfer functions....
void PhotonPairSelector::FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ) {

  pt = -999.;
  decayZ = -999.;

  // loop over all GEN particles and look for status 1 photons
  for(UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle* p = fMCParticles->At(i);
    if( !(p->Is(MCParticle::kH)) ) continue;
    pt=p->Pt();
    decayZ = p->DecayVertex().Z();
    break;
  }
  
  return;
 }

// this routine looks for the idx of the run-range
Int_t PhotonPairSelector::FindRunRangeIdx(UInt_t run) {
  Int_t runRange=-1;
  for(UInt_t iRun = 0; iRun<fRunStart.size(); ++iRun) {
    if( run >= fRunStart[iRun] && run <= fRunEnd[iRun]) {
      runRange = (Int_t) iRun;
      return runRange;
    }
  }
  return runRange;
}


Double_t PhotonPairSelector::GetDataEnCorr(Int_t runRange, PhotonTools::CiCBaseLineCats cat) {
  switch( cat ) {
  case PhotonTools::kCiCCat1:
    return fDataEnCorr_EB_hR9[runRange];
  case PhotonTools::kCiCCat2:
    return fDataEnCorr_EB_lR9[runRange];
  case PhotonTools::kCiCCat3:
    return fDataEnCorr_EE_hR9[runRange];
  case PhotonTools::kCiCCat4:
    return fDataEnCorr_EE_lR9[runRange];
  default:
    return 1.;
  }
}


Double_t PhotonPairSelector::GetMCSmearFac(PhotonTools::CiCBaseLineCats cat) {
  switch( cat ) {
  case PhotonTools::kCiCCat1:
    return fMCSmear_EB_hR9;
  case PhotonTools::kCiCCat2:
    return fMCSmear_EB_lR9;
  case PhotonTools::kCiCCat3:
    return fMCSmear_EE_hR9;
  case PhotonTools::kCiCCat4:
    return fMCSmear_EE_lR9;
  default:
    return 1.;
  }
}

Float_t PhotonPairSelector::GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2) {
  
  bool ph1IsEB = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat2);
  bool ph2IsEB = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat2);

  bool ph1IsHR9 = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat3);
  bool ph2IsHR9 = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat3);
  
  if( ph1IsEB && ph2IsEB )
    return ( ph1IsHR9 && ph2IsHR9 ? 0. : 1.);
  
  return ( ph1IsHR9 && ph2IsHR9 ? 2. : 3.);
}

void PhotonPairSelectorPhoton::SetVars(const Photon *p, const MCParticle *m) {
      e = p->E();
      pt = p->Pt();
      eta = p->Eta();
      phi = p->Phi();
      r9 = p->R9();
      e3x3 = p->E33();
      e5x5 = p->E55();
      const SuperCluster *s = p->SCluster();
      sce = s->Energy();
      scrawe = s->RawEnergy();
      scpse = s->PreshowerEnergy();
      sceta = s->Eta();
      scphi = s->Phi();
      scnclusters = s->ClusterSize();
      scnhits = s->NHits();
      hovere = p->HadOverEm();
      sigietaieta = p->CoviEtaiEta();
      isbarrel = (s->AbsEta()<1.5);
      isr9reco = (isbarrel && r9>0.94) || (!isbarrel && r9>0.95);
      isr9cat = (r9>0.94);
      
      if (isbarrel) {
        if (isr9cat) phcat = 1;
        else phcat = 2;
      }
      else {
        if (isr9cat) phcat = 3;
        else phcat = 4;
      }
      
      const BasicCluster *b = s->Seed();
      sigiphiphi = TMath::Sqrt(b->CoviPhiiPhi());
      if (isnan(sigiphiphi)) sigiphiphi = -99.;
      covietaiphi = b->CoviEtaiPhi();
      if (isnan(covietaiphi)) covietaiphi = -99.;
      emax = b->EMax();
      e2nd = b->E2nd();
      etop = b->ETop();
      ebottom = b->EBottom();
      eleft = b->ELeft();
      eright = b->ERight();
      e1x3 = b->E1x3();
      e3x1 = b->E3x1();
      e1x5 = b->E1x5();
      e2x2 = b->E2x2();
      e4x4 = b->E4x4();
      e2x5max = b->E2x5Max();
      e2x5top = b->E2x5Top();
      e2x5bottom = b->E2x5Bottom();
      e2x5left = b->E2x5Left();
      e2x5right = b->E2x5Right();

      //initialize photon energy corrections if needed
      if (!PhotonFix::initialised()) {
        PhotonFix::initialise("4_2",std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFix.dat")).Data()));  
      }
      
      PhotonFix fix(e,sceta,scphi,r9);
      const Float_t dval = -99.;
      ecor = fix.fixedEnergy();
      ecorerr = fix.sigmaEnergy();
      if (isbarrel) {
       etac = fix.etaC();
       etas = fix.etaS();
       etam = fix.etaM();
       phic = fix.phiC();
       phis = fix.phiS();
       phim = fix.phiM();
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
       xz = fix.xZ();
       xc = fix.xC();
       xs = fix.xS();
       xm = fix.xM();
       yz = fix.yZ();
       yc = fix.yC();
       ys = fix.yS();
       ym = fix.yM();
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
      }
      else {
        ispromptgen = kFALSE;
        gene = -99.;
        genpt = -99.;
        geneta = -99.;
        genphi = -99.;
      }
      
}