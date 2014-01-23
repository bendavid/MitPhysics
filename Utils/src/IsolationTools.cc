// $Id: IsolationTools.cc,v 1.37 2013/02/23 14:53:07 mingyang Exp $

#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

ClassImp(mithep::IsolationTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::TrackIsolation(const Track *p, Double_t extRadius, Double_t intRadius,
                                        Double_t ptLow, Double_t maxVtxZDist, 
                                        const Collection<Track> *tracks) 
{
  //Computes the Track Isolation: Summed Transverse Momentum of all tracks inside an 
  //annulus around the electron seed track.  

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<tracks->GetEntries();i++) {   
    Double_t tmpPt = tracks->At(i)->Pt();
    Double_t deltaZ = fabs(p->Z0() - tracks->At(i)->Z0());

    //ignore the track if it is below the pt threshold
    if (tmpPt < ptLow) 
      continue;    
    //ingore the track if it is too far away in Z
    if (deltaZ > maxVtxZDist) 
      continue;
           
    Double_t dr = MathUtils::DeltaR(p->Phi(),p->Eta(),tracks->At(i)->Phi(), tracks->At(i)->Eta());
    //add the track pt if it is inside the annulus
    if ( dr < extRadius && 
	 dr >= intRadius ) {
      ptSum += tmpPt;
    }
  }
  return ptSum;  
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::EcalIsolation(const SuperCluster *sc, Double_t coneSize, Double_t etLow, 
                                       const Collection<BasicCluster> *basicClusters) 
{
  //Computes the Ecal Isolation: Summed Transverse Energy of all Basic Clusters inside a  
  //cone around the electron, excluding those that are inside the electron super cluster.

  Double_t ecalIsol=0.;
  const BasicCluster *basicCluster= 0;
  for (UInt_t i=0; i<basicClusters->GetEntries();i++) {    
    basicCluster = basicClusters->At(i);    
    Double_t basicClusterEnergy    = basicCluster->Energy();
    Double_t basicClusterEta  = basicCluster->Eta();
    Double_t basicClusterEt   = basicClusterEnergy*sin(2*atan(exp(basicClusterEta)));           

    if (basicClusterEt > etLow) {            
      bool inSuperCluster = false;	  
      
      // loop over the basic clusters of the supercluster
      // to make sure that the basic cluster is not inside
      // the super cluster. We exclude those.
      for (UInt_t j=0; j<sc->ClusterSize(); j++) {
        const BasicCluster *tempBasicClusterInSuperCluster = sc->Cluster(j);	
        if (tempBasicClusterInSuperCluster == basicCluster) {
          inSuperCluster = true;	    
        }
      }
      
      if (!inSuperCluster) {	    
        Double_t dr = MathUtils::DeltaR(sc->Phi(), sc->Eta(),
                                      basicCluster->Phi(),basicCluster->Eta());
        if(dr < coneSize) {
          ecalIsol += basicClusterEt;
        }
      }
    }
  } 
  return ecalIsol;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::CaloTowerHadIsolation(const ThreeVector *p, Double_t extRadius, 
                                               Double_t intRadius, Double_t etLow, 
                                               const Collection<CaloTower> *caloTowers) 
{
  //Computes the CaloTower Had Et Isolation: Summed Hadronic Transverse Energy of all Calo Towers 
  //inside an annulus around the electron super cluster position.

  Double_t sumEt = 0;
  for (UInt_t i=0; i<caloTowers->GetEntries();i++) {    
    Double_t caloTowerEt = caloTowers->At(i)->HadEt();
    Double_t dr = MathUtils::DeltaR(caloTowers->At(i)->Phi(), caloTowers->At(i)->Eta(),
                                  p->Phi(), p->Eta());
    if (dr < extRadius && dr > intRadius && caloTowerEt > etLow) {
      sumEt += caloTowerEt;
    }
  }
  return sumEt;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::CaloTowerEmIsolation(const ThreeVector *p, Double_t extRadius, 
                                              Double_t intRadius, Double_t etLow, 
                                              const Collection<CaloTower> *caloTowers) 
{
  //Computes the CaloTower Em Et Isolation: Summed Hadronic Transverse Energy of all Calo Towers 
  //inside an annulus around the electron super cluster position.

  Double_t sumEt = 0;
  for (UInt_t i=0; i<caloTowers->GetEntries();i++) {    
    Double_t caloTowerEt = caloTowers->At(i)->EmEt();
    Double_t dr = MathUtils::DeltaR(caloTowers->At(i)->Phi(), caloTowers->At(i)->Eta(),
                                   p->Phi(), p->Eta());
    if (dr < extRadius && dr > intRadius && caloTowerEt > etLow) {
      sumEt += caloTowerEt;
    }
  }
  return sumEt;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFRadialMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, 
                                               Double_t ptMin, Double_t extRadius)
{
  //Computes the PF Radial Muon Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track with a dR weighting

  double RadialIso = 0;
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);

    // exclude muon
    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    
    // inner cone veto for tracks
    if (dr < 0.01) continue;

    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    // exclude PFMuons and PFElectrons within the cone
    if (pf->HasTrk() &&
       (pf->PFType() == PFCandidate::eElectron || 
        pf->PFType() == PFCandidate::eMuon)) continue;

    // add the pf pt if it is inside the extRadius 
    if (dr < extRadius) RadialIso += pf->Pt() * (1.0 - 3.0*dr);
  }
  return RadialIso;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, 
                                      	 const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
				     	 Double_t extRadius, Double_t intRadiusGamma, Double_t intRadius)
{
  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);

    // exclude muon
    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    
    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->HasTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }
      
    // inner cone veto for gammas
    if (pf->PFType() == PFCandidate::eGamma && dr < intRadiusGamma) continue;
    
    // inner cone veto for tracks
    if (dr < intRadius) continue;

    // add the pf pt if it is inside the extRadius 
    if (dr < extRadius) ptSum += pf->Pt();
  
  }
  return ptSum;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands,
                                      	 const MuonCol *goodMuons, const ElectronCol *goodElectrons, 
                                      	 const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
				      	 Double_t extRadius, Double_t intRadiusGamma, Double_t intRadius)
{
  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);

    // exclude muon
    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    
    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->HasTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }
      
    // inner cone veto for gammas
    if (pf->PFType() == PFCandidate::eGamma && dr < intRadiusGamma) continue;
    
    // inner cone veto for tracks
    if (dr < intRadius) continue;

    // add the pf pt if it is inside the extRadius and outside the intRadius
    if (dr < extRadius ) {
      Bool_t isLepton = kFALSE;
      if(goodMuons){
        for (UInt_t nl=0; nl<goodMuons->GetEntries();nl++) {
	  const Muon *m = goodMuons->At(nl);
          if(pf->TrackerTrk() && m->TrackerTrk() &&
	     pf->TrackerTrk() == m->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }
      if(goodElectrons && isLepton == kFALSE){
        for (UInt_t nl=0; nl<goodElectrons->GetEntries();nl++) {
	  const Electron *e = goodElectrons->At(nl);
          if(pf->TrackerTrk() && e->TrackerTrk() &&
	     pf->TrackerTrk() == e->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
          if(pf->GsfTrk() && e->GsfTrk() &&
	     pf->GsfTrk() == e->GsfTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }
      if (isLepton == kTRUE) continue;
      ptSum += pf->Pt();            
    }  
  }
  return ptSum;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                                      	     const Vertex *vertex, Double_t delta_z, Double_t ptMin,
				     	     Double_t extRadius, Double_t intRadius, Int_t PFCandidateType) 
{

  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);
    
    //select only specified PFCandidate types
    if (PFCandidateType >= 0) {
      if (pf->PFType() != PFCandidateType) continue;
    }

    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    if(pf->GsfTrk() && p->GsfTrk() &&
       pf->GsfTrk() == p->GsfTrk()) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->BestTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }

           
    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    // add the pf pt if it is inside the extRadius and outside the intRadius
    if ( dr < extRadius && 
	 dr >= intRadius ) {

      //EtaStrip Veto for Gamma 
      if (pf->PFType() == PFCandidate::eGamma && fabs(p->Eta() - pf->Eta()) < 0.025) continue;

      //InnerCone (One Tower = dR < 0.07) Veto for non-gamma neutrals
      if (!pf->HasTrk() && pf->PFType() == PFCandidate::eNeutralHadron
          && MathUtils::DeltaR(p->Mom(), pf->Mom()) < 0.07 ) continue; 

      ptSum += pf->Pt();            

    }
  }
  return ptSum;
}


//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                      	     		     const MuonCol *goodMuons, const ElectronCol *goodElectrons,
					     const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
			     		     Double_t extRadius, Double_t intRadius, Int_t PFCandidateType)
{

  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);
    
    //select only specified PFCandidate types
    if (PFCandidateType >= 0) {
      if (pf->PFType() != PFCandidateType) continue;
    }

    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    if(pf->GsfTrk() && p->GsfTrk() &&
       pf->GsfTrk() == p->GsfTrk()) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->BestTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }
           
    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    // add the pf pt if it is inside the extRadius and outside the intRadius
    if ( dr < extRadius && 
	 dr >= intRadius ) {

      //EtaStrip Veto for Gamma 
      if (pf->PFType() == PFCandidate::eGamma && fabs(p->Eta() - pf->Eta()) < 0.025) continue;

      //InnerCone (One Tower = dR < 0.07) Veto for non-gamma neutrals
      if (!pf->HasTrk() && pf->PFType() == PFCandidate::eNeutralHadron
          && MathUtils::DeltaR(p->Mom(), pf->Mom()) < 0.07 ) continue; 

      Bool_t isLepton = kFALSE;
      if(goodMuons){
        for (UInt_t nl=0; nl<goodMuons->GetEntries();nl++) {
	  const Muon *m = goodMuons->At(nl);
          if(pf->TrackerTrk() && m->TrackerTrk() &&
	     pf->TrackerTrk() == m->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }
      if(goodElectrons && isLepton == kFALSE){
        for (UInt_t nl=0; nl<goodElectrons->GetEntries();nl++) {
	  const Electron *e = goodElectrons->At(nl);
          if(pf->TrackerTrk() && e->TrackerTrk() &&
	     pf->TrackerTrk() == e->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
          if(pf->GsfTrk() && e->GsfTrk() &&
	     pf->GsfTrk() == e->GsfTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }

      if (isLepton == kTRUE) continue;
      ptSum += pf->Pt();            

    }
  }
  return ptSum;
}
//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFElectronIsolation2012(const Electron *ele, const Vertex *vertex, 
                                                 const PFCandidateCol *PFCands, 
                                                 const PileupEnergyDensityCol *PileupEnergyDensity,
                                                 ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaTarget,
                                                 const ElectronCol *goodElectrons,
                                                 const MuonCol *goodMuons, Double_t dRMax, Bool_t isDebug){

  Double_t tmpChargedIso_DR       = 0;
  Double_t tmpGammaIso_DR         = 0;
  Double_t tmpNeutralHadronIso_DR = 0;

  for (UInt_t p=0; p<PFCands->GetEntries();p++) {   
    const PFCandidate *pf = PFCands->At(p);

    //exclude the electron itself
    //if(pf->GsfTrk() && ele->GsfTrk() &&
    //   pf->GsfTrk() == ele->GsfTrk()) continue;
    //if(pf->TrackerTrk() && ele->TrackerTrk() &&
    //   pf->TrackerTrk() == ele->TrackerTrk()) continue;      

    //************************************************************
    // New Isolation Calculations
    //************************************************************
    Double_t dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());

    if (dr < 1.0) {
      Bool_t IsLeptonFootprint = kFALSE;
      //************************************************************
      // Lepton Footprint Removal
      //************************************************************            
      for (UInt_t q=0; q < goodElectrons->GetEntries() ; ++q) {
	//if pf candidate matches an electron passing ID cuts, then veto it
	if(pf->GsfTrk() && goodElectrons->At(q)->GsfTrk() &&
	   pf->GsfTrk() == goodElectrons->At(q)->GsfTrk()) IsLeptonFootprint = kTRUE;
	if(pf->TrackerTrk() && goodElectrons->At(q)->TrackerTrk() &&
	   pf->TrackerTrk() == goodElectrons->At(q)->TrackerTrk()) IsLeptonFootprint = kTRUE;
	//if pf candidate lies in veto regions of electron passing ID cuts, then veto it
	if(pf->BestTrk() && fabs(goodElectrons->At(q)->SCluster()->Eta()) >= 1.479 
           && MathUtils::DeltaR(goodElectrons->At(q)->Mom(), pf->Mom()) < 0.015) IsLeptonFootprint = kTRUE;
	if(pf->PFType() == PFCandidate::eGamma && fabs(goodElectrons->At(q)->SCluster()->Eta()) >= 1.479 &&
	   MathUtils::DeltaR(goodElectrons->At(q)->Mom(), pf->Mom()) < 0.08) IsLeptonFootprint = kTRUE;
      }
      for (UInt_t q=0; q < goodMuons->GetEntries() ; ++q) {
	//if pf candidate matches an muon passing ID cuts, then veto it
	if(pf->TrackerTrk() && goodMuons->At(q)->TrackerTrk() &&
	   pf->TrackerTrk() == goodMuons->At(q)->TrackerTrk()) IsLeptonFootprint = kTRUE;
	//if pf candidate lies in veto regions of muon passing ID cuts, then veto it
	if(pf->BestTrk() && MathUtils::DeltaR(goodMuons->At(q)->Mom(), pf->Mom()) < 0.01) IsLeptonFootprint = kTRUE;
      }

      if (!IsLeptonFootprint) {
	Bool_t passVeto = kTRUE;
	//Charged
	 if(pf->BestTrk()) {	  	   
	   // CMS DOESN"T WANT THIS
	   //if (!(fabs(pf->BestTrk()->DzCorrected(*vertex) - ele->BestTrk()->DzCorrected(*vertex)) < 0.2)) passVeto = kFALSE;
	   //************************************************************
	   // Veto any PFmuon, or PFEle
	   if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = kFALSE;
	   //************************************************************
	   //************************************************************
	   // Footprint Veto
	   if (fabs(ele->SCluster()->Eta()) >= 1.479 && dr < 0.015) passVeto = kFALSE;
	   //************************************************************
	   if (passVeto) {
	     if (dr < dRMax) tmpChargedIso_DR += pf->Pt();
	   } //pass veto	  
	 }
	 //Gamma
	 else if (pf->PFType() == PFCandidate::eGamma) {
	   //************************************************************
	   // Footprint Veto
	   if (fabs(ele->SCluster()->Eta()) >= 1.479) {
             if (dr < 0.08) passVeto = kFALSE;
	   }
	   //************************************************************
	   
	   if (passVeto) {
	     if (dr < dRMax) tmpGammaIso_DR += pf->Pt();
	   }
	 }
	 //NeutralHadron
	 else {
           if (dr < dRMax) tmpNeutralHadronIso_DR += pf->Pt();
	 }
      } //not lepton footprint
    } //in 1.0 dr cone
  } //loop over PF candidates

  Double_t Rho = 0;
  if (!(TMath::IsNaN(PileupEnergyDensity->At(0)->Rho()) || isinf(PileupEnergyDensity->At(0)->Rho()))) Rho = PileupEnergyDensity->At(0)->Rho();
  
  Double_t IsoVar_ChargedIso_DR = tmpChargedIso_DR/ele->Pt();
  Double_t IsoVar_NeutralIso_DR = tmpGammaIso_DR + tmpNeutralHadronIso_DR;
  //   Careful here, we have kEleNeutralIso04 only for now            ****added kEleGammaIso03 and kEleNeutralHadronIso03 --heng june 16th 2012
  if(dRMax != 0.4) assert(0);
  double EA = ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralIso04, ele->SCluster()->Eta(), EffectiveAreaTarget);

  IsoVar_NeutralIso_DR = TMath::Max((IsoVar_NeutralIso_DR - Rho*EA)/ele->Pt(), 0.0);
  
  if(isDebug == kTRUE){
    printf("Iso(ch, em, nh), EA, RelIso = (%f, %f, %f), %f, %f\n",tmpChargedIso_DR,tmpGammaIso_DR,tmpNeutralHadronIso_DR,
           EA,IsoVar_ChargedIso_DR + IsoVar_NeutralIso_DR);
  }

  return (IsoVar_ChargedIso_DR + IsoVar_NeutralIso_DR);
}

//--------------------------------------------------------------------------------------------------

Double_t IsolationTools::PFElectronIsolation2012LepTag(const Electron *ele, const Vertex *vertex, 
						       const PFCandidateCol *PFCands, 
						       const PileupEnergyDensityCol *PileupEnergyDensity,
						       ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaTarget,
						       const ElectronCol *goodElectrons,
						       const MuonCol *goodMuons, Double_t dRMax, Bool_t isDebug){




  Double_t tmpChargedIso_DR       = 0;
  Double_t tmpGammaIso_DR         = 0;
  Double_t tmpNeutralHadronIso_DR = 0;

  for (UInt_t p=0; p<PFCands->GetEntries();p++) {   
    const PFCandidate *pf = PFCands->At(p);

    //************************************************************
    // New Isolation Calculations
    //************************************************************
    Double_t dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());
    if (dr < 1.0) {
      Bool_t IsLeptonFootprint = kFALSE;
      //************************************************************
      // Lepton Footprint Removal
      //************************************************************   
      if(pf->PFType() == PFCandidate::eGamma){
	for (UInt_t i=0; i<PFCands->GetEntries();i++) { 
	  const PFCandidate *pf2 = PFCands->At(i);
	  if(pf2->PFType()==PFCandidate::eElectron && pf->SCluster() == pf2->SCluster())IsLeptonFootprint = kTRUE;
	}
      }
         
      for (UInt_t q=0; q < goodElectrons->GetEntries() ; ++q) {
	//if pf candidate matches an electron passing ID cuts, then veto it
	if(pf->GsfTrk() && goodElectrons->At(q)->GsfTrk() &&
	   pf->GsfTrk() == goodElectrons->At(q)->GsfTrk()) IsLeptonFootprint = kTRUE;
	if(pf->TrackerTrk() && goodElectrons->At(q)->TrackerTrk() &&
	   pf->TrackerTrk() == goodElectrons->At(q)->TrackerTrk()) IsLeptonFootprint = kTRUE;
	//if pf candidate lies in veto regions of electron passing ID cuts, then veto it
	if(pf->BestTrk() && fabs(goodElectrons->At(q)->SCluster()->Eta()) >= 1.479 
           && MathUtils::DeltaR(goodElectrons->At(q)->Mom(), pf->Mom()) < 0.015) IsLeptonFootprint = kTRUE;
	if(pf->PFType() == PFCandidate::eGamma && fabs(goodElectrons->At(q)->SCluster()->Eta()) >= 1.479 &&
	   MathUtils::DeltaR(goodElectrons->At(q)->Mom(), pf->Mom()) < 0.08) IsLeptonFootprint = kTRUE;
      }
      for (UInt_t q=0; q < goodMuons->GetEntries() ; ++q) {
	//if pf candidate matches an muon passing ID cuts, then veto it
	if(pf->TrackerTrk() && goodMuons->At(q)->TrackerTrk() &&
	   pf->TrackerTrk() == goodMuons->At(q)->TrackerTrk()) IsLeptonFootprint = kTRUE;
	//if pf candidate lies in veto regions of muon passing ID cuts, then veto it
	if(pf->BestTrk() && MathUtils::DeltaR(goodMuons->At(q)->Mom(), pf->Mom()) < 0.01) IsLeptonFootprint = kTRUE;
      }

      if (!IsLeptonFootprint) {
	Bool_t passVeto = kTRUE;
	//Charged
	 if(pf->BestTrk()) {	  	   
	   // CMS DOESN"T WANT THIS
	   //if (!(fabs(pf->BestTrk()->DzCorrected(*vertex) - ele->BestTrk()->DzCorrected(*vertex)) < 0.2)) passVeto = kFALSE;
	   //************************************************************
	   // Veto any PFmuon, or PFEle
	   if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = kFALSE;
	   //************************************************************
	   //************************************************************
	   // Footprint Veto
	   if (fabs(ele->SCluster()->Eta()) >= 1.479 && dr < 0.015) passVeto = kFALSE;
	   //************************************************************
	   if (passVeto) {
	     if (dr < dRMax) {
	       if ( false ) printf("Charged PFCand Pt = %f    eta= %f       phi= %f     , dr= %f    charge=%f   \n",pf->Pt(),pf->Eta(),pf->Phi(),dr, pf->Charge());
	       tmpChargedIso_DR += pf->Pt();
	     } //pass dr
	   } //pass veto	  
	 }
	 //Gamma
	 else if (pf->PFType() == PFCandidate::eGamma) {
	   //************************************************************
	   // Footprint Veto
	   if (fabs(ele->SCluster()->Eta()) >= 1.479) {
             if (dr < 0.08) passVeto = kFALSE;
	   }
	   //************************************************************
	   
	   if (passVeto) {
	     if(pf->HasSCluster()){
	       if ( false ) printf("Gamma PFCand Pt = %f    eta= %f       phi= %f     , dr= %f   \n",pf->Pt(),pf->SCluster()->Eta(),pf->Phi(),dr);//ming
	     }
	     //if (dr < dRMax ) tmpGammaIso_DR += pf->Pt();
	     //if (dr < dRMax && ((pf->HasSCluster() && pf->SCluster()!=ele->SCluster()) || !pf->HasSCluster())) tmpGammaIso_DR += pf->Pt();//ming
	     if(dr < dRMax && !(ele->GsfTrk()->NExpectedHitsInner()>0 && pf->MvaGamma() > 0.99 && pf->HasSCluster() && ele->SCluster() == pf->SCluster())) tmpGammaIso_DR += pf->Pt();//ming:HZZ
	   }
	 }
	 //NeutralHadron
	 else {
	   if ( false ) printf("NeutralHadron PFCand Pt = %f    eta= %f       phi= %f     , dr= %f   \n",pf->Pt(),pf->Eta(),pf->Phi(),dr);
           if (dr < dRMax) tmpNeutralHadronIso_DR += pf->Pt();
	 }
      } //not lepton footprint
    } //in 1.0 dr cone
  } //loop over PF candidates

  Double_t Rho = 0;
  assert(PileupEnergyDensity);
  //if (!(TMath::IsNaN(PileupEnergyDensity->At(0)->Rho()) || isinf(PileupEnergyDensity->At(0)->Rho()))) Rho = PileupEnergyDensity->At(0)->Rho();//ming
  Rho = PileupEnergyDensity->At(0)->RhoKt6PFJets();//ming

  Double_t IsoVar_ChargedIso_DR = tmpChargedIso_DR/ele->Pt();
  Double_t IsoVar_NeutralIso_DR = tmpGammaIso_DR + tmpNeutralHadronIso_DR;
  // Careful here, we have kEleNeutralIso04 only for now            ****added kEleGammaIso03 and kEleNeutralHadronIso03 --heng june 16th 2012
  double EA = ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso03, ele->SCluster()->Eta(), EffectiveAreaTarget) + ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso03, ele->SCluster()->Eta(), EffectiveAreaTarget) ;
  IsoVar_NeutralIso_DR = TMath::Max((IsoVar_NeutralIso_DR - Rho*EA)/ele->Pt(), 0.0);
  
  if(isDebug == kTRUE){
    printf("Iso(ch, em, nh), EA, RelIso = (%f, %f, %f), %f, %f\n",tmpChargedIso_DR,tmpGammaIso_DR,tmpNeutralHadronIso_DR,
           EA,IsoVar_ChargedIso_DR + IsoVar_NeutralIso_DR);
  }

  return (IsoVar_ChargedIso_DR + IsoVar_NeutralIso_DR);
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::BetaM(const TrackCol *tracks, const Muon *p, const Vertex *vertex, 
                               Double_t ptMin, Double_t  delta_z, Double_t extRadius,
			       Double_t intRadius){

  if(!tracks) return 1.0;
  if(tracks->GetEntries() <= 0) return 1.0;

  double Pt_jets_X = 0. ;
  double Pt_jets_Y = 0. ;
  double Pt_jets_X_tot = 0. ;
  double Pt_jets_Y_tot = 0. ;

  for(int i=0;i<int(tracks->GetEntries());i++){
    const Track *pTrack = tracks->At(i);

    if(pTrack && p->TrackerTrk() &&
       pTrack == p->TrackerTrk()) continue;

    if(pTrack->Pt() <= ptMin) continue;

    Double_t dr = MathUtils::DeltaR(pTrack->Mom(),p->Mom());
    if ( dr < extRadius && dr >= intRadius ) {
      Pt_jets_X_tot += pTrack->Px();
      Pt_jets_Y_tot += pTrack->Py();  
      double pDz = TMath::Abs(pTrack->DzCorrected(*vertex));
      if(pDz < delta_z){
        Pt_jets_X += pTrack->Px();
        Pt_jets_Y += pTrack->Py();
      }
    }
  }

  if(sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot) > 0)
    return sqrt(Pt_jets_X*Pt_jets_X + Pt_jets_Y*Pt_jets_Y) / sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot);

  return 1.0;
}



//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::BetaMwithPUCorrection(const PFCandidateCol *PFNoPileUp, const PFCandidateCol *PFPileUp, const Muon *p, Double_t extRadius){

  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t ptSumChHadnoPU =0.;  
  Double_t etsumNeuHad =0;
  Double_t ptsumGamma = 0;
  Double_t ptsumPU= 0;
  Double_t ptSum= 0;
  for (UInt_t i=0; i<PFNoPileUp->GetEntries();i++) { 
    const PFCandidate *pf = PFNoPileUp->At(i);
    // exclude muon
    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    //exclude PFCands outside of cone 0.4
    Double_t dr1 = MathUtils::DeltaR(p->Mom(), pf->Mom());
    if (dr1 > extRadius) continue;
    
    if (pf->PFType() == PFCandidate::eHadron) ptSumChHadnoPU += pf->Pt();
    
    //inner veto for neutrals
    if  (dr1<0.01) continue;
    
    //pt threshold for neutrals;
    if (pf->Pt()<0.5) continue;
    
    if (pf->PFType() == PFCandidate::eNeutralHadron) etsumNeuHad += pf->Pt();
    else if (pf->PFType() == PFCandidate::eGamma) ptsumGamma += pf->Pt(); 
  }

  for (UInt_t i=0; i<PFPileUp->GetEntries();i++) { 
    const PFCandidate *pf = PFPileUp->At(i);
    // exclude muon
    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;
    //exclude PFCands outside of cone 0.4
    Double_t dr2 = MathUtils::DeltaR(p->Mom(), pf->Mom());
    if (dr2 > extRadius) continue;
    
    //inner veto for neutrals
    if  (dr2<0.01) continue;
    
    //pt threshold for neutrals;
    if (pf->Pt()<0.5) continue;    
    
    if (pf->PFType() == PFCandidate::eHadron || pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) ptsumPU += pf->Pt();
  }
  
  //   for (UInt_t i=0; i<PFPUCands->GetEntries(); ++i){
  //     const PFCandidate *pfpu = PFPUCand->At(i);
  //     ptsumPU += pfpu->Pt();
  //   }
  //printf("ChHad    Gamma     Neu      PU    \n");
  // printf("%f       %f       %f        %f     \n",ptSumChHadnoPU,ptsumGamma,etsumNeuHad,ptsumPU);
  ptSum=ptSumChHadnoPU+ TMath::Max(0.,ptsumGamma+etsumNeuHad - 0.5*ptsumPU);
  //printf("ming muon check isolation is %f \n",ptSum);
  
  return ptSum;
}






//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::BetaE(const TrackCol *tracks, const Electron *p, const Vertex *vertex, 
                               Double_t ptMin, Double_t  delta_z, Double_t extRadius,
			       Double_t intRadius){

  if(!tracks) return 1.0;
  if(tracks->GetEntries() <= 0) return 1.0;

  double Pt_jets_X = 0. ;
  double Pt_jets_Y = 0. ;
  double Pt_jets_X_tot = 0. ;
  double Pt_jets_Y_tot = 0. ;

  for(int i=0;i<int(tracks->GetEntries());i++){
    const Track *pTrack = tracks->At(i);

    if(pTrack && p->TrackerTrk() &&
       pTrack == p->TrackerTrk()) continue;

    if(pTrack && p->GsfTrk() &&
       pTrack == p->GsfTrk()) continue;

    if(pTrack->Pt() <= ptMin) continue;

    Double_t dr = MathUtils::DeltaR(pTrack->Mom(),p->Mom());
    if ( dr < extRadius && dr >= intRadius ) {
      Pt_jets_X_tot += pTrack->Px();
      Pt_jets_Y_tot += pTrack->Py();  
      double pDz = TMath::Abs(pTrack->DzCorrected(*vertex));
      if(pDz < delta_z){
        Pt_jets_X += pTrack->Px();
        Pt_jets_Y += pTrack->Py();
      }
    }
  }

  if(sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot) > 0)
    return sqrt(Pt_jets_X*Pt_jets_X + Pt_jets_Y*Pt_jets_Y) / sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot);

  return 1.0;
}


// method added by F.Stoeckli: computes the track isolation with NO constrint on the OV-track compatibility
Double_t IsolationTools::TrackIsolationNoPV(const mithep::Particle* p, const BaseVertex* bsp, 
					    Double_t extRadius, 
					    Double_t intRadius, 
					    Double_t ptLow, 
					    Double_t etaStrip,
					    Double_t maxD0,
					    mithep::TrackQuality::EQuality quality,
					    const mithep::Collection<mithep::Track> *tracks,
                                            UInt_t maxNExpectedHitsInner,
                                            const mithep::DecayParticleCol *conversions) {
  
  // loop over all tracks
  Double_t tPt = 0.;
  //std::cout<<"      *** TrackIso:"<<std::endl;
  for(UInt_t i=0; i<tracks->GetEntries(); ++i) {
    const Track* t = tracks->At(i);
    if(t->Pt()>1. && false) {
      std::cout<<"     "<<i<<"    pt = "<<t->Pt()<<"  ("<<ptLow<<")"<<std::endl;
      std::cout<<"                d0 = "<<fabs(t->D0Corrected( *bsp) )<<"  ("<<maxD0<<")"<<std::endl;
      //std::cout<<"              conv ? "<<PhotonTools::MatchedConversion(t,conversions,bsp)<<std::endl;
      std::cout<<"                dR = "<<MathUtils::DeltaR(t->Mom(),p->Mom())<<"  ("<<extRadius<<","<<intRadius<<")"<<std::endl;
      std::cout<<"             dEta  = "<<fabs(t->Eta()-p->Eta())<<"   ("<<etaStrip<<")"<<std::endl;
    }
    if ( t->Pt() < ptLow ) continue;
    if ( ! t->Quality().Quality(quality) ) continue;
    // only check for beamspot if available, otherwise ignore cut
    if ( bsp && fabs(t->D0Corrected( *bsp) ) > maxD0) continue;
    if (t->NExpectedHitsInner()>maxNExpectedHitsInner) continue;
    if (conversions && PhotonTools::MatchedConversion(t,conversions,bsp)) continue;
    Double_t dR   = MathUtils::DeltaR(t->Mom(),p->Mom());
    Double_t dEta = fabs(t->Eta()-p->Eta());
    if(dR < extRadius && dR > intRadius && dEta > etaStrip) tPt += t->Pt();
  }
  return tPt;
}


Double_t IsolationTools::CiCTrackIsolation(const mithep::Photon* p, 
					   const BaseVertex* theVtx, 
					   Double_t extRadius, 
					   Double_t intRadius, 
					   Double_t ptLow, 
					   Double_t etaStrip,
					   Double_t maxD0,
					   Double_t maxDZ,
					   const mithep::Collection<mithep::Track> *tracks,
					   unsigned int* worstVtxIndex,
					   const mithep::Collection<mithep::Vertex> *vtxs,
					   const mithep::Collection<mithep::Electron> *eles,
					   bool print,
					   double* ptmax, double* dRmax) {
  
  UInt_t numVtx = 1;
  const BaseVertex* iVtx = theVtx;
  if( vtxs ) { 
    numVtx = vtxs->GetEntries();
    if (numVtx > 0)
      iVtx = vtxs->At(0);
    else
      return 0.;
  }
  
  // NEW for Electron T&P: need to remove the electron Gsf Track (applied if eles != NULL)
  const Track* theGsfTrack = NULL;
  if ( eles ) {
    // find the electron that matches the Photon SC
    for(UInt_t j=0; j<eles->GetEntries(); ++j) {
      if ( eles->At(j)->SCluster() == p->SCluster() ) {
	if( eles->At(j)->HasTrackerTrk() )
	  theGsfTrack = eles->At(j)->TrackerTrk();
	break;
      }
    }
  }

  if(print) {
    std::cout<<" Testing photon with"<<std::endl;
    std::cout<<"             Et  = "<<p->Et()<<std::endl;
    std::cout<<"             Eta = "<<p->Eta()<<std::endl;
    std::cout<<"             Phi = "<<p->Phi()<<std::endl;
  }

  Double_t iIso = 0.;
  Double_t maxIso = 0.;
  
  if(worstVtxIndex)
    *worstVtxIndex=0;

  double t_ptmax = 0.;
  double t_dRmax = 0.;

  for(UInt_t i=0; i<numVtx; ++i) {

    if(i>0) iVtx = vtxs->At(i);


    if(print) {
      std::cout<<"   Vertex #"<<i<<std::endl;
      std::cout<<"       with X = "<<iVtx->X()<<std::endl;
      std::cout<<"       with Y = "<<iVtx->Y()<<std::endl;
      std::cout<<"       with Z = "<<iVtx->Z()<<std::endl;
    }

    Photon* phTemp = new Photon(*p);

    // RESET CALO_POS!
    phTemp->SetCaloPosXYZ(p->SCluster()->Point().X(),p->SCluster()->Point().Y(),p->SCluster()->Point().Z());

    // compute the ph momentum with respect to this Vtx
    //FourVectorM phMom = p->MomVtx(iVtx->Position());
    FourVectorM phMom = phTemp->MomVtx(iVtx->Position());

    delete phTemp;

    if(print) {
      std::cout<<"         photon has changed to:"<<std::endl;
      std::cout<<"             Et  = "<<phMom.Et()<<std::endl;
      std::cout<<"             eta = "<<phMom.Eta()<<std::endl;
      std::cout<<"             Phi = "<<phMom.Phi()<<std::endl;
    }
    
    iIso = 0.;
    for(UInt_t j=0; j<tracks->GetEntries(); ++j) {
      const Track* t = tracks->At(j);
      if( theGsfTrack && t == theGsfTrack ) continue;

      //Double_t dR   = MathUtils::DeltaR(t->Mom(),p->Mom());
      //Double_t dEta = TMath::Abs(t->Eta()-p->Eta());

      Double_t dR   = MathUtils::DeltaR(t->Mom(),phMom);
      Double_t dEta = TMath::Abs(t->Eta()-phMom.Eta());

      if(print && t->Pt()>1. && false) {
	  std::cout<<"              passing track #"<<j<<std::endl;
	  std::cout<<"                          pt = "<<t->Pt()<<std::endl;
	  std::cout<<"                         eta = "<<t->Eta()<<std::endl;
	  std::cout<<"                         phi = "<<t->Phi()<<std::endl;
	  std::cout<<"                          d0 = "<<fabs(t->D0Corrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dZ = "<<fabs(t->DzCorrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dR = "<<dR<<std::endl;
	  std::cout<<"                        dEta = "<<dEta<<std::endl;
	  std::cout<<"                          vx = "<<t->X0()<<std::endl;
	  std::cout<<"                          vy = "<<t->Y0()<<std::endl;
	  std::cout<<"                          vz = "<<t->Z0()<<std::endl;
      }

      if ( t->Pt() < ptLow ) continue;
      // only check for beamspot if available, otherwise ignore cut
      if ( fabs(t->D0Corrected( *iVtx )) > maxD0) continue;
      if ( fabs(t->DzCorrected( *iVtx )) > maxDZ) continue;


      if(dR < extRadius && dR > intRadius && dEta > etaStrip) {
	iIso += t->Pt();

	if(t->Pt() > t_ptmax) {
	  t_ptmax=t->Pt();
	  t_dRmax=dR;
	}

	if(print && t->Pt()>1.) {
	  std::cout<<"              passing track #"<<j<<std::endl;
	  std::cout<<"                          pt = "<<t->Pt()<<std::endl;
	  std::cout<<"                         eta = "<<t->Eta()<<std::endl;
	  std::cout<<"                         phi = "<<t->Phi()<<std::endl;
	  std::cout<<"                          d0 = "<<fabs(t->D0Corrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dZ = "<<fabs(t->DzCorrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dR = "<<dR<<std::endl;
	  std::cout<<"                        dEta = "<<dEta<<std::endl;
	  std::cout<<"                          vx = "<<t->X0()<<std::endl;
	  std::cout<<"                          vy = "<<t->Y0()<<std::endl;
	  std::cout<<"                          vz = "<<t->Z0()<<std::endl;
	  std::cout<<"                             new  tIso = "<<iIso<<std::endl;
	}
      }
    }
    if ( iIso > maxIso ) {
      maxIso = iIso;
      if(worstVtxIndex)
	*worstVtxIndex=i;
    }
  }

  if(ptmax) (*ptmax)=t_ptmax;
  if(dRmax) (*dRmax)=t_dRmax;

  if(print) {
    if(worstVtxIndex)
      std::cout<<"   max TrkIso is given by Vtx #"<<*worstVtxIndex<<" with an amount of tIso = "<<maxIso<<std::endl;
    else
      std::cout<<"   max TrkIso is given by Vtx #0 with an amount of tIso = "<<maxIso<<std::endl;
  } 
  return maxIso;
}

//ChargedIso_selvtx_DR0To0p001=IsolationTools::PFChargedIsolation(p, SelVtx, 0.01, 0, 0.0, 0.0, 0.1, 0.2,fPFCands);

Double_t IsolationTools::PFChargedIsolation(const mithep::Photon *p, 
                                         const BaseVertex *theVtx, 
                                         Double_t extRadius,
                                         Double_t intRadius,
                                         const PFCandidateCol *PFCands,
                                         unsigned int* worstVtxIndex,
                                         const mithep::Collection<mithep::Vertex> *vtxs,
                                         bool print)
{
  
  UInt_t numVtx = 1;
  
  const BaseVertex* iVtx = theVtx;
  
  if( vtxs ) { 
    numVtx = vtxs->GetEntries();
    if (numVtx > 0)
      iVtx = vtxs->At(0);
    else
      return 0.;
  }
    
  if(print) {
    std::cout<<" Testing photon with"<<std::endl;
    std::cout<<"             Et  = "<<p->Et()<<std::endl;
    std::cout<<"             Eta = "<<p->Eta()<<std::endl;
    std::cout<<"             Phi = "<<p->Phi()<<std::endl;
  }
  
  Double_t iIso = 0.;
  Double_t maxIso = 0.;
  
  if(worstVtxIndex)
    *worstVtxIndex=0;
  
  for(UInt_t i=0; i<numVtx; ++i) {
    
    if(i>0) iVtx = vtxs->At(i);
    
    
    if(print) {
      std::cout<<"   Vertex #"<<i<<std::endl;
      std::cout<<"       with X = "<<iVtx->X()<<std::endl;
      std::cout<<"       with Y = "<<iVtx->Y()<<std::endl;
      std::cout<<"       with Z = "<<iVtx->Z()<<std::endl;
    }
        
    ThreeVector photondir = ThreeVector(p->SCluster()->Point()) - iVtx->Position();
        
    iIso = 0.;
    
    for(UInt_t j=0; j<PFCands->GetEntries(); ++j) {
      const PFCandidate *pf= PFCands->At(j);

      if( pf->HasTrackerTrk() && pf->PFType()==PFCandidate::eHadron ) {
	const Track* t = pf->TrackerTrk();
        
	Double_t dR   = MathUtils::DeltaR(*pf,photondir);
	//Double_t dEta = TMath::Abs(pf->Eta()-photondir.Eta());
	
        if (dR<0.02) continue;
        if (dR<intRadius) continue;
        if (dR>extRadius) continue;
        
	if(print && pf->Pt()>1. && false) {
	  std::cout<<"              passing track #"<<j<<std::endl;
	  std::cout<<"                          pt = "<<pf->Pt()<<std::endl;
	  std::cout<<"                         eta = "<<pf->Eta()<<std::endl;
	  std::cout<<"                         phi = "<<pf->Phi()<<std::endl;
	  std::cout<<"                          d0 = "<<fabs(t->D0Corrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dZ = "<<fabs(t->DzCorrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dR = "<<dR<<std::endl;
	  //std::cout<<"                        dEta = "<<dEta<<std::endl;
	  std::cout<<"                          vx = "<<t->X0()<<std::endl;
	  std::cout<<"                          vy = "<<t->Y0()<<std::endl;
	  std::cout<<"                          vz = "<<t->Z0()<<std::endl;
	}
	
	
	// only check for beamspot if available, otherwise ignore cut
	if ( fabs(t->D0Corrected( *iVtx )) > 0.1) continue;
	if ( fabs(t->DzCorrected( *iVtx )) > 0.2) continue;
	
	iIso += pf->Pt();
	
      }
    }
    
    if ( iIso > maxIso ) {
      maxIso = iIso;
      if(worstVtxIndex)
	*worstVtxIndex=i;
    }
  }
  

  if(print) {
    if(worstVtxIndex)
      std::cout<<"   max TrkIso is given by Vtx #"<<*worstVtxIndex<<" with an amount of tIso = "<<maxIso<<std::endl;
    else
      std::cout<<"   max TrkIso is given by Vtx #0 with an amount of tIso = "<<maxIso<<std::endl;
  }
  return maxIso;
}

//
Float_t IsolationTools::PFChargedCount(const mithep::Photon* p, 
				   const BaseVertex* theVtx, 
				   Double_t extRadius, 
				   Double_t intRadius, 
				   Double_t ptLow, 
				   Double_t etaStrip,
				   Double_t maxD0,
				   Double_t maxDZ,
				   const PFCandidateCol *PFCands,
				   unsigned int* worstVtxIndex,
				   const mithep::Collection<mithep::Vertex> *vtxs,
				   const mithep::Collection<mithep::Electron> *eles,
				   bool print,
				   double* ptmax,
				   double* dRmax) {
  
  UInt_t numVtx = 1;
  
  const BaseVertex* iVtx = theVtx;
  
  if( vtxs ) { 
    numVtx = vtxs->GetEntries();
    if (numVtx > 0)
      iVtx = vtxs->At(0);
    else
      return 0.;
  }
  
  // NEW for Electron T&P: need to remove the electron Gsf Track (applied if eles != NULL)
  const Track* theGsfTrack = NULL;
  if ( eles ) {
    // find the electron that matches the Photon SC
    for(UInt_t j=0; j<eles->GetEntries(); ++j) {
      if ( eles->At(j)->SCluster() == p->SCluster() ) {
	if( eles->At(j)->HasTrackerTrk() )
	  theGsfTrack = eles->At(j)->TrackerTrk();
	break;
      }
    }
  }
  
  if(print) {
    std::cout<<" Testing photon with"<<std::endl;
    std::cout<<"             Et  = "<<p->Et()<<std::endl;
    std::cout<<"             Eta = "<<p->Eta()<<std::endl;
    std::cout<<"             Phi = "<<p->Phi()<<std::endl;
  }
  
  Double_t iIso = 0.;
  Double_t maxIso = 0.;

  Float_t iNumParticles = 0.;
  Float_t maxNumParticles = 0.;
  
  if(worstVtxIndex)
    *worstVtxIndex=0;
  
  double t_ptmax = 0.;
  double t_dRmax = 0.;
  
  for(UInt_t i=0; i<numVtx; ++i) {
    
    if(i>0) iVtx = vtxs->At(i);
    
    
    if(print) {
      std::cout<<"   Vertex #"<<i<<std::endl;
      std::cout<<"       with X = "<<iVtx->X()<<std::endl;
      std::cout<<"       with Y = "<<iVtx->Y()<<std::endl;
      std::cout<<"       with Z = "<<iVtx->Z()<<std::endl;
    }
    
    Photon* phTemp = new Photon(*p);
    
    // RESET CALO_POS! //ming: why?
    phTemp->SetCaloPosXYZ(p->SCluster()->Point().X(),p->SCluster()->Point().Y(),p->SCluster()->Point().Z());
    
    // compute the ph momentum with respect to this Vtx
    FourVectorM phMom = phTemp->MomVtx(iVtx->Position());
    
    delete phTemp;
    
    if(print) {
      std::cout<<"         photon has changed to:"<<std::endl;
      std::cout<<"             Et  = "<<phMom.Et()<<std::endl;
      std::cout<<"             eta = "<<phMom.Eta()<<std::endl;
      std::cout<<"             Phi = "<<phMom.Phi()<<std::endl;
    }
    
    iIso = 0.;
    iNumParticles = 0.;
    
    for(UInt_t j=0; j<PFCands->GetEntries(); ++j) {
      const PFCandidate *pf= PFCands->At(j);
      if(pf->HasTrk() && (pf->PFType()==PFCandidate::eHadron || pf->PFType()==PFCandidate::eElectron || pf->PFType()==PFCandidate::eMuon)){
	const Track* t = pf->BestTrk();
        if(pf->PFType()==PFCandidate::eElectron && pf->HasGsfTrk()){t = pf->GsfTrk();}
	if(!(pf->PFType()==PFCandidate::eElectron) && pf->HasTrackerTrk()){t = pf->TrackerTrk();}
	
	if( theGsfTrack && t == theGsfTrack ) continue;
	
	Double_t dR   = MathUtils::DeltaR(pf->Mom(),phMom);
	Double_t dEta = TMath::Abs(pf->Eta()-phMom.Eta());
	
	if(print && pf->Pt()>1. && false) {
	  std::cout<<"              passing track #"<<j<<std::endl;
	  std::cout<<"                          pt = "<<pf->Pt()<<std::endl;
	  std::cout<<"                         eta = "<<pf->Eta()<<std::endl;
	  std::cout<<"                         phi = "<<pf->Phi()<<std::endl;
	  std::cout<<"                          d0 = "<<fabs(t->D0Corrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dZ = "<<fabs(t->DzCorrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dR = "<<dR<<std::endl;
	  std::cout<<"                        dEta = "<<dEta<<std::endl;
	  std::cout<<"                          vx = "<<t->X0()<<std::endl;
	  std::cout<<"                          vy = "<<t->Y0()<<std::endl;
	  std::cout<<"                          vz = "<<t->Z0()<<std::endl;
	}
	
	if ( pf->Pt() < ptLow ) continue;
	
	// only check for beamspot if available, otherwise ignore cut
	if ( fabs(t->D0Corrected( *iVtx )) > maxD0) continue;
	if ( fabs(t->DzCorrected( *iVtx )) > maxDZ) continue;
	
	
	if(dR < extRadius && dR > intRadius && dEta > etaStrip) {
	  iIso += pf->Pt();
	  iNumParticles += 1;

	  if(pf->Pt() > t_ptmax) {
	    t_ptmax=pf->Pt();
	    t_dRmax=dR;
	  }
	  
	  if(print && pf->Pt()>1.) {
	    std::cout<<"              passing track #"<<j<<std::endl;
	    std::cout<<"                          pt = "<<pf->Pt()<<std::endl;
	    std::cout<<"                         eta = "<<pf->Eta()<<std::endl;
	    std::cout<<"                         phi = "<<pf->Phi()<<std::endl;
	    std::cout<<"                          d0 = "<<fabs(t->D0Corrected( *iVtx ))<<std::endl;
	    std::cout<<"                          dZ = "<<fabs(t->DzCorrected( *iVtx ))<<std::endl;
	    std::cout<<"                          dR = "<<dR<<std::endl;
	    std::cout<<"                        dEta = "<<dEta<<std::endl;
	    std::cout<<"                          vx = "<<t->X0()<<std::endl;
	    std::cout<<"                          vy = "<<t->Y0()<<std::endl;
	    std::cout<<"                          vz = "<<t->Z0()<<std::endl;
	    std::cout<<"                             new  tIso = "<<iIso<<std::endl;
	  }
	}
      }
    }
    
    if ( iIso > maxIso ) {
      maxIso = iIso;
      maxNumParticles = iNumParticles;
      
      if(worstVtxIndex)
	*worstVtxIndex=i;
    }
  }
  
  if(ptmax) (*ptmax)=t_ptmax;
  if(dRmax) (*dRmax)=t_dRmax;

  if(print) {
    if(worstVtxIndex)
      std::cout<<"   max TrkIso is given by Vtx #"<<*worstVtxIndex<<" with an amount of tIso = "<<maxIso<<std::endl;
    else
      std::cout<<"   max TrkIso is given by Vtx #0 with an amount of tIso = "<<maxIso<<std::endl;
  }
  return maxNumParticles;
}

Double_t IsolationTools::PFGammaIsolation(const mithep::Photon *p, 
                                         Double_t extRadius,
                                         Double_t intRadius,
                                         const PFCandidateCol *PFCands)
{
 
  Double_t iso = 0.;
  
  ThreeVector photondir;
  //Bool_t setdir = kFALSE;
  
  for (UInt_t ipfc = 0; ipfc<PFCands->GetEntries(); ++ipfc) {
    const PFCandidate *pfc = PFCands->At(ipfc);
    if (pfc->PFType()!=PFCandidate::eGamma) continue;

    ThreeVector photondir = ThreeVector(p->SCluster()->Point() - pfc->SourceVertex());
    
    Double_t dR   = MathUtils::DeltaR(*pfc,photondir);
    Double_t dEta = TMath::Abs(pfc->Eta()-photondir.Eta());
    
    Bool_t isbarrel = p->SCluster()->AbsEta()<1.5;
    
    if (isbarrel && dEta<0.015) continue;
    if (!isbarrel && dR<0.07) continue;
    if (dR<intRadius) continue;
    if (dR>extRadius) continue;
    if (pfc->HasSCluster() && pfc->SCluster()==p->SCluster()) continue;
    
    iso += pfc->Pt();
    
  }
  
  return iso;
  
  
}

Double_t IsolationTools::PFNeutralHadronIsolation(const mithep::Photon *p, 
                                         Double_t extRadius,
                                         Double_t intRadius,
                                         const PFCandidateCol *PFCands)
{
 
  Double_t iso = 0.;
  
  ThreeVector photondir;
  Bool_t setdir = kFALSE;
  
  for (UInt_t ipfc = 0; ipfc<PFCands->GetEntries(); ++ipfc) {
    const PFCandidate *pfc = PFCands->At(ipfc);
    if (pfc->PFType()!=PFCandidate::eNeutralHadron) continue;
    if (!setdir) {
      photondir = ThreeVector(p->SCluster()->Point() - pfc->SourceVertex());
      setdir = kTRUE;
    }
    
    Double_t dR   = MathUtils::DeltaR(*pfc,photondir);

    if (dR<intRadius) continue;
    if (dR>extRadius) continue;
    
    iso += pfc->Pt();
    
  }
  
  return iso;
  
  
}
