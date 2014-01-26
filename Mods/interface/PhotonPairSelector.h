//--------------------------------------------------------------------------------------------------
// M.Yang 2011/10/12
// $Id: PhotonPairSelector.h,v 1.37 2013/12/09 17:55:50 bendavid Exp $
//
// PhotonPairSelector
//
// Authors: F. Stoeckli
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONPAIRSELECTOR_H
#define MITPHYSICS_MODS_PHOTONPAIRSELECTOR_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"

#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/EGEnergyCorrector.h"

#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"

#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"

class TNtuple;
class TRandom3;
class TH1D;

namespace mithep 
{  
  class PhotonPairSelector : public BaseMod
  {
  public:
    PhotonPairSelector(const char *name  = "PhotonPairSelector", 
		       const char *title = "Selecting PhotonPairs");
    
    ~PhotonPairSelector();

    enum PhotonSelection {
      kNoPhSelection = 0,
      kCiCPhSelection,
      kCiCPFPhSelection,   
      kCiCPFPhSelectionNoPFChargedIso,   
      kMVAPhSelection, //MVA
      kMITPhSelection,
      kMITPFPhSelection,
      kMITPFPhSelectionNoEcal,
      kMITPFPhSelectionNoEcalNoPFChargedIso
    };
    enum VertexSelection {
      kStdVtxSelection = 0,
      kCiCVtxSelection,
      kMITVtxSelection,
      kCiCMVAVtxSelection,
      kMetSigVtxSelection
    };

    class EnergyCorrection {
      
      public:
        
        EnergyCorrection() : fCat(PhotonTools::kEBlowEtaGold), fMinRun(0), fMaxRun(0), fMinEt(-99.), fMaxEt(-99.), fCorr(1.) {}
        EnergyCorrection(PhotonTools::eScaleCats cat, UInt_t minrun, UInt_t maxrun, Double_t minet, Double_t maxet, Double_t corr) : 
          fCat(cat), fMinRun(minrun), fMaxRun(maxrun), fMinEt(minet), fMaxEt(maxet), fCorr(corr) {}
          
        PhotonTools::eScaleCats fCat;
        UInt_t fMinRun;
        UInt_t fMaxRun;
        Double_t fMinEt;
        Double_t fMaxEt;
        Double_t fCorr;
      
    };
    
    
    // outsourced to MVATools (IdMVAType)  (fab)
/*     enum IdMVA { */
/*       k2011IdMVA = 0, */
/*       k2012IdMVA_globe */
/*     }; */
   

    // outsourced to PhotonTools (fab)
/*     enum ShowerShape { */
/*       k2011ShowerShape  = 0, */
/*       k2012ShowerShape */
/*     }; */

    // setting all the input Names
    void                SetInputPhotonsName(const char *n){ fPhotonBranchName= n;        }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetTrackName(const char *n)       { fTrackBranchName = n;        }
    void                SetElectronName(const char *n)    { fElectronName = n;           }
    void                SetConversionName(const char *n)  { fConversionName = n;         }
    void                SetPUDensityName(const char *n)   { fPileUpDenName = n;          }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetMCParticle(const char *n)      { fMCParticleName = n;         }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }
    void                SetPFCandName(const char *n)      { fPFCandName = n;             }

    // set the type of selection
    void                SetPhotonSelType(const char *type){ fPhotonSelType    = type;    }
    void                SetVertexSelType(const char *type){ fVertexSelType    = type;    }
    void                SetIdMVAType(const char *type)    { fIdMVAType        = type;    } 
    void                SetShowerShapeType(const char *type) { fShowerShapeType        = type;    } 

    // get/set the Names for the output Photon Collection
    const char         *GetOutputName()             const { return fGoodPhotonsName;     }   
    void                SetOutputName(const char *n)      { fGoodPhotonsName=n;          }    

    const char         *GetOutputVtxName()          const { return fChosenVtxName;       }
    void                SetOutputVtxName(const char* n)   { fChosenVtxName = n;          }


    // set basic Cut variables (FOR PRE-SELECTION)
    void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
    void                SetAbsEtaMax(Double_t eta)        { fPhotonEtaMax    = eta;	 }
    
    void                SetLeadingPtMin(Double_t pt)      { fLeadingPtMin = pt;          }
    void                SetTrailingPtMin(Double_t pt)     { fTrailingPtMin = pt;         }

    // is Data Or Not?
    void                SetIsData (Bool_t b) { fIsData = b;}

    
    void AddEnCorrFromFile(TString filename);
    
    // methods to set the MC smearing/energy correction values
    void                AddEnCorrPerRun( PhotonTools::eScaleCats cat, UInt_t sRun, UInt_t eRun, Double_t corr) {

      fDataEnCorr.push_back(EnergyCorrection(cat,sRun,eRun,-99.,-99.,corr));
    };

    void                AddEnCorrPerRunEtDep( PhotonTools::eScaleCats cat, UInt_t sRun, UInt_t eRun, Double_t minet, Double_t maxet, Double_t corr) {

      fDataEnCorr.push_back(EnergyCorrection(cat,sRun,eRun,minet,maxet,corr));
    };    
    
    void                SetMCSmearFactors(Double_t _EBlowEta_hR9,
					  Double_t _EBlowEta_lR9,
                                          Double_t _EBhighEta_hR9, 
                                          Double_t _EBhighEta_lR9,
					  Double_t _EElowEta_hR9, 
					  Double_t _EElowEta_lR9,
                                          Double_t _EEhighEta_hR9,
                                          Double_t _EEhighEta_lR9) {
      fMCSmear_EBlowEta_hR9 = _EBlowEta_hR9;
      fMCSmear_EBlowEta_lR9 = _EBlowEta_lR9;
      fMCSmear_EBhighEta_hR9 = _EBhighEta_hR9;
      fMCSmear_EBhighEta_lR9 = _EBhighEta_lR9;      
      fMCSmear_EElowEta_hR9 = _EElowEta_hR9;
      fMCSmear_EElowEta_lR9 = _EElowEta_lR9;
      fMCSmear_EEhighEta_hR9 = _EEhighEta_hR9;
      fMCSmear_EEhighEta_lR9 = _EEhighEta_lR9;      
    };

    // special routine in ordetr to use different smearing to compute mass-error for di-photon MVA input
    void                SetMCSmearFactorsMVA(Double_t _EBlowEta_hR9,
					     Double_t _EBlowEta_lR9,
					     Double_t _EBhighEta_hR9, 
					     Double_t _EBhighEta_lR9,
					     Double_t _EElowEta_hR9, 
					     Double_t _EElowEta_lR9,
					     Double_t _EEhighEta_hR9,
					     Double_t _EEhighEta_lR9) {

      fMCSmearMVA_EBlowEta_hR9 = _EBlowEta_hR9;
      fMCSmearMVA_EBlowEta_lR9 = _EBlowEta_lR9;
      fMCSmearMVA_EBhighEta_hR9 = _EBhighEta_hR9;
      fMCSmearMVA_EBhighEta_lR9 = _EBhighEta_lR9;      
      fMCSmearMVA_EElowEta_hR9 = _EElowEta_hR9;
      fMCSmearMVA_EElowEta_lR9 = _EElowEta_lR9;
      fMCSmearMVA_EEhighEta_hR9 = _EEhighEta_hR9;
      fMCSmearMVA_EEhighEta_lR9 = _EEhighEta_lR9;      
    };

    void SetMCStochasticPivot(Double_t pivot_lowEta_hR9,
                              Double_t pivot_lowEta_lR9,
                              Double_t pivot_highEta_hR9,
                              Double_t pivot_highEta_lR9) {
     
      fMCStochasticPivot_EBlowEta_hR9 = pivot_lowEta_hR9;
      fMCStochasticPivot_EBlowEta_lR9 = pivot_lowEta_lR9;
      fMCStochasticPivot_EBhighEta_hR9 = pivot_highEta_hR9;
      fMCStochasticPivot_EBhighEta_lR9 = pivot_highEta_lR9;
      
    }
    
    void SetMCStochasticRho(Double_t rho_lowEta_hR9,
                              Double_t rho_lowEta_lR9,
                              Double_t rho_highEta_hR9,
                              Double_t rho_highEta_lR9) {
     
      fMCStochasticRho_EBlowEta_hR9 = rho_lowEta_hR9;
      fMCStochasticRho_EBlowEta_lR9 = rho_lowEta_lR9;
      fMCStochasticRho_EBhighEta_hR9 = rho_highEta_hR9;
      fMCStochasticRho_EBhighEta_lR9 = rho_highEta_lR9;
      
    }    

    void SetMCStochasticPhi(Double_t phi_lowEta_hR9,
                              Double_t phi_lowEta_lR9,
                              Double_t phi_highEta_hR9,
                              Double_t phi_highEta_lR9) {
     
      fMCStochasticPhi_EBlowEta_hR9 = phi_lowEta_hR9;
      fMCStochasticPhi_EBlowEta_lR9 = phi_lowEta_lR9;
      fMCStochasticPhi_EBhighEta_hR9 = phi_highEta_hR9;
      fMCStochasticPhi_EBhighEta_lR9 = phi_highEta_lR9;
      
    }    
    

    void                SetApplyEleVeto(bool a)            { fApplyEleVeto  = a; }
    void                SetInvertElectronVeto(Bool_t b)   { fInvertElectronVeto = b;     }          
    void                DoDataEneCorr(bool a)           { fDoDataEneCorr = a; }
    void                DoMCSmear(bool a)               { fDoMCSmear     = a; }
    void                DoMCEneSmear(bool a)            { fDoMCEneSmear = a; }
    void                DoEneErrSmear(bool a)           { fDoEneErrSmear = a; }
    void                UseSpecialSmearForDPMVA(bool a) { fUseSpecialSmearForDPMVA     = a; }
    void                SetStochasticSmear(bool a)      { fStochasticSmear = a; }
    
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name) { fGoodElectronName = name; }
    void                SetUseSingleLegConversions(Bool_t b) { fUseSingleLegConversions = b; }
    void                SetDoRegression(Bool_t b)         { fDoRegression = b; }
    void                SetEtaCorrections(const TH1D *h)  { fEtaCorrections = h; }
    void                SetBdtCutBarrel(Float_t x)        { fbdtCutBarrel = x; }
    void                SetBdtCutEndcap(Float_t x)        { fbdtCutEndcap = x; }
 
    void                SetDoShowerShapeScaling(Bool_t b)  { fDoShowerShapeScaling = b; }
    
    void                SetJetsName(const char *n)       { fJetsName = n;              }    
        
    void                SetRhoType(RhoUtilities::RhoType type) { fRhoType = type ; }

    void                SetApplyLeptonTag(bool a)            {  fApplyLeptonTag = a; }

    void                SetLeptonTagElectronsName(TString name) { fLeptonTagElectronsName = name; }
    void                SetLeptonTagMuonsName    (TString name) { fLeptonTagMuonsName     = name; }

    void                SetR9Rescale(bool dor9rescale, double p0b, double p1b, double  p0e,  double  p1e) {fdor9rescale = dor9rescale; fp0b = p0b; fp1b = p1b; fp0e = p0e; fp1e = p1e;}

  protected:
    void                Process();
    void                SlaveBegin();

    // private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z, Float_t& mass);
    Double_t            GetDataEnCorr(UInt_t run, const Photon *p);
    Double_t            GetMCSmearFac(PhotonTools::eScaleCats cat,                    bool useSpecialSmear = false);   // last flag in case of special smearing for error computation
    
    Double_t            GetMCSmearFacStochastic(Photon *p) const;
    
    Float_t             GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2);

    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fElectronName;
    TString             fGoodElectronName;
    TString             fConversionName;
    TString             fPFConversionName;    
    TString             fTrackBranchName;
    TString             fPileUpDenName;    
    TString             fPVName;
    TString             fBeamspotName;
    TString             fPFCandName;
    TString             fMCParticleName;
    TString             fPileUpName;
    TString             fJetsName;
    TString             fPFMetName;
    
    TString             fGoodPhotonsName;      //name of exported "good photon" collection
    TString             fChosenVtxName;        //name of exported "chosen Vtx"  collection

    TString             fLeptonTagElectronsName;
    TString             fLeptonTagMuonsName;
    
    // Selection Types
    TString             fPhotonSelType;
    TString             fVertexSelType;
    PhotonSelection     fPhSelType;
    VertexSelection     fVtxSelType;    

    // Id Type
    TString             fIdMVAType;
    MVATools::IdMVAType fIdType;
    //IdMVA               fIdType;

    // showershape
    TString                                fShowerShapeType;
    PhotonTools::ShowerShapeScales         fSSType;

    // Basic Pre-Selection kinematics
    Double_t            fPhotonPtMin;          // min pt cut fro PRE-SELECTION!
    Double_t            fPhotonEtaMax;         // max eta cut for PRE-SELECTION!
    
    Double_t            fLeadingPtMin;
    Double_t            fTrailingPtMin;

    // is it Data or MC?
    Bool_t              fIsData;
    
    // in case there's some PV pre-selection
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;
    Bool_t              fUseSingleLegConversions;

    Bool_t              fStochasticSmear;

    const PhotonCol              *fPhotons;
    const ElectronCol            *fElectrons;
    const ElectronCol            *fGoodElectrons;    
    const DecayParticleCol       *fConversions;
    const DecayParticleCol       *fPFConversions;
    const TrackCol               *fTracks;
    const PileupEnergyDensityCol *fPileUpDen;
    const VertexCol              *fPV;
    const BeamSpotCol            *fBeamspot;
    const PFCandidateCol         *fPFCands;
    const MCParticleCol          *fMCParticles;
    const PileupInfoCol          *fPileUp;    
    const JetCol                 *fJets;
    const PFMetCol               *fPFMet;
    const ElectronCol            *fLeptonTagElectrons;
    const MuonCol                *fLeptonTagMuons;
       
    //vector of energy scale corrections
    //format of tuple is:
    //category, minrun, maxrun, minet, maxet, correction
    std::vector<EnergyCorrection> fDataEnCorr;
    
    Double_t              fMCSmear_EBlowEta_hR9;
    Double_t              fMCSmear_EBlowEta_lR9;
    Double_t              fMCSmear_EBhighEta_hR9;
    Double_t              fMCSmear_EBhighEta_lR9;    
    Double_t              fMCSmear_EElowEta_hR9;
    Double_t              fMCSmear_EElowEta_lR9;
    Double_t              fMCSmear_EEhighEta_hR9;
    Double_t              fMCSmear_EEhighEta_lR9;    

    //  special Smear factors for usage for diphoton MVA input, incase differrent from std smearing
    Double_t              fMCSmearMVA_EBlowEta_hR9;
    Double_t              fMCSmearMVA_EBlowEta_lR9;
    Double_t              fMCSmearMVA_EBhighEta_hR9;
    Double_t              fMCSmearMVA_EBhighEta_lR9;    
    Double_t              fMCSmearMVA_EElowEta_hR9;
    Double_t              fMCSmearMVA_EElowEta_lR9;
    Double_t              fMCSmearMVA_EEhighEta_hR9;
    Double_t              fMCSmearMVA_EEhighEta_lR9;    
    
    Double_t              fMCStochasticPivot_EBlowEta_hR9;
    Double_t              fMCStochasticPivot_EBlowEta_lR9;
    Double_t              fMCStochasticPivot_EBhighEta_hR9;
    Double_t              fMCStochasticPivot_EBhighEta_lR9;     
    
    Double_t              fMCStochasticRho_EBlowEta_hR9;
    Double_t              fMCStochasticRho_EBlowEta_lR9;
    Double_t              fMCStochasticRho_EBhighEta_hR9;
    Double_t              fMCStochasticRho_EBhighEta_lR9; 
    
    Double_t              fMCStochasticPhi_EBlowEta_hR9;
    Double_t              fMCStochasticPhi_EBlowEta_lR9;
    Double_t              fMCStochasticPhi_EBhighEta_hR9;
    Double_t              fMCStochasticPhi_EBhighEta_lR9;     


    // pointer to RNG ionstance for smearing
    TRandom3*             fRng;
    EGEnergyCorrector     fEgCor;
    Bool_t                fDoRegression;
    TString               fPhFixString;
    TString               fPhFixFile;
    TString               fRegWeights; 
  		         
    const TH1D           *fEtaCorrections;
    
    // --------------------------------
    // some streagin flags, not adjustable yet (FIX-ME)
    bool                  fDoDataEneCorr;
    bool                  fDoMCSmear;
    bool                  fDoMCEneSmear;
    bool                  fDoEneErrSmear;
    bool                  fUseSpecialSmearForDPMVA;   // if set to true, the special smearing numbers set in fMCSmearMVA_* are used to compute the mass-errors (input to diphoton MVA)
    bool                  fDoVtxSelection;
    bool                  fApplyEleVeto;
    Bool_t                fInvertElectronVeto; //=true (invert ele veto, for cic sel only atm)
    
    // --------------------------------
    bool                  fApplyLeptonTag;   

    //MVA
    int                   fVariableType_2011;
    TString               fEndcapWeights_2011;
    TString               fBarrelWeights_2011;
    int                   fVariableType_2012_globe;
    TString               fEndcapWeights_2012_globe;
    TString               fBarrelWeights_2012_globe;  
    MVATools              fTool;
    Float_t               fbdtCutBarrel;
    Float_t               fbdtCutEndcap;

    VertexTools           fVtxTools;
    
    Bool_t                fDoShowerShapeScaling; 
       	                  
    Bool_t                fDoMCErrScaling;
    Double_t              fMCErrScaleEB;
    Double_t              fMCErrScaleEE;    
    UInt_t                fRegressionVersion;
    	                  
    Bool_t                fRelativePtCuts;
    
    MVAMet                fMVAMet;

    //R9 rescale
    bool                  fdor9rescale;                   
    double                fp0b;                          
    double                fp1b;                          
    double                fp0e;                          
    double                fp1e;                          
    
    RhoUtilities::RhoType fRhoType;

    ClassDef(PhotonPairSelector, 1) // Photon identification module
  };
}
#endif
