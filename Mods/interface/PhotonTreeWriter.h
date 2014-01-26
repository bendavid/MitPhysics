//--------------------------------------------------------------------------------------------------
// $Id: PhotonTreeWriter.h,v 1.51 2013/12/12 17:14:08 veverka Exp $
//
// PhotonTreeWriter
//
// Authors: J. Bendavid, J. Veverka
//--------------------------------------------------------------------------------------------------

/**
 * TODO
 * - Factor out taggers (VHLep, VHHad, ttH) into separate classes
 * - Factor out transient diphoton-final-state data (phHard, phSoft, selvtx, ..)
 *   into a separate class.
 */


#ifndef MITPHYSICS_MODS_PHOTONTREEWRITER_H
#define MITPHYSICS_MODS_PHOTONTREEWRITER_H


#include <utility>  // provides std::pair
#include <vector>
#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/Particle.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/MVAVBF.h"

#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

class TNtuple;
class TRandom3;

namespace mithep 
{

  class PhotonTreeWriterDiphotonEvent;
  
  class PhotonTreeWriterVtx
  {
    
    public:
      void SetVars(const Vertex *v, const Photon *p1, const Photon *p2, const PFCandidateCol *pfcands, Int_t idx, Int_t numvtx, Float_t genvtxz);
      
    private:
      Int_t n;
      Int_t nvtx;
      Float_t ptgg;
      Float_t phigg;
      Float_t etagg;
      Float_t pxgg;
      Float_t pygg;
      Float_t pzgg;
      Float_t mgg;
      Float_t zgen;
      Float_t x;
      Float_t y;
      Float_t z;
      Float_t pt;
      Float_t phi;
      Float_t eta;
      Float_t px;
      Float_t py;
      Float_t pz;
      Float_t sumpt;
      Float_t sumptsq;
      Int_t nchalltoward;
      Int_t nchalltransverse;
      Int_t nchallaway;
      Int_t nchcuttoward;
      Int_t nchcuttransverse;
      Int_t nchcutaway;
  };
  
  template <int NClus=16>//ming: what's template
  class PhotonTreeWriterPhoton
  {
    public:  
    void SetVars(const Photon *p,
                 const DecayParticle *c,
                 const Electron *ele,
                 const SuperCluster *pfsc,
                 const MCParticle *m,
                 PhotonFix &phfixph,
                 PhotonFix &phfixele,
                 const TrackCol* trackCol,
                 const VertexCol* vtxCol,
                 const PFCandidateCol* candCol,
                 Double_t _tRho,
                 Bool_t fillclusterarrays,
                 const PhotonCol* ps,
                 const SuperClusterCol* scs,
                 const ElectronCol* els,
                 const DecayParticleCol *convs,
                 const BaseVertex *bs,bool isdata,
		 bool dor9rescale, double p0b, double p1b, double  p0e,  double  p1e,
                 Bool_t applyElectronVeto=kTRUE,
                 const Vertex* realVtx = NULL);
      Float_t Ecor()    const { return ecor;    };
      Float_t Ecorerr() const { return ecorerr; };
      Float_t Ecorele()    const { return ecorele;    };
      Float_t Ecoreleerr() const { return ecoreleerr; };
      Float_t Pt() const {return pt;};
      
    private:  
      UChar_t hasphoton;
      UInt_t  index;
      Float_t e;
      Float_t pt;
      Float_t eta;
      Float_t phi;
      Float_t r9;
      Float_t e3x3;
      Float_t e5x5;
      UInt_t  scindex;
      Float_t sce;
      Float_t scrawe;
      Float_t scpse;
      Float_t scpssigmaxx;
      Float_t scpssigmayy;
      Float_t sceta;
      Float_t scphi;
      UInt_t  scnclusters;
      UInt_t  scnhits;
      Float_t scetawidth;
      Float_t scphiwidth;
      Float_t hovere;
      Float_t hoveretower;
      Float_t sigietaieta;
      Bool_t  isbarrel;
      Bool_t  isr9reco;
      Bool_t  isr9cat;
      Char_t  phcat;
      Float_t eerr;
      Float_t eerrsmeared;
      Float_t esmearing;
      Float_t eerrsmearing;
      Float_t epresmearscale;
      Float_t escale;
      Float_t stochasticrho;
      Float_t stochasticphi;
      Float_t stochasticpivot;
      Float_t idmva;
      Float_t ecalisodr03;
      Float_t hcalisodr03;
      Float_t trkisohollowdr03;
      Float_t ecalisodr04;
      Float_t hcalisodr04;
      Float_t trkisohollowdr04;
      Float_t trackiso1;
      Float_t trackiso2;
      Float_t combiso1;
      Float_t combiso2;
      Bool_t passeleveto;
      
      //quantities from seed basic cluster      
      Float_t eseed;
      Float_t etaseed;
      Float_t phiseed;
      Int_t   ietaseed;
      Int_t   iphiseed;
      Int_t   ixseed;
      Int_t   iyseed;
      Float_t etacryseed;
      Float_t phicryseed;
      Float_t xcryseed;
      Float_t ycryseed;
      Float_t thetaaxisseed;
      Float_t phiaxisseed;
      Float_t sigietaietaseed;
      Float_t sigiphiphiseed;
      Float_t covietaiphiseed;
      Float_t e3x3seed;
      Float_t e5x5seed;      
      Float_t emaxseed;
      Float_t e2ndseed;
      Float_t etopseed;
      Float_t ebottomseed;
      Float_t eleftseed;
      Float_t erightseed;      
      Float_t e1x3seed;
      Float_t e3x1seed;
      Float_t e1x5seed;
      Float_t e2x2seed;
      Float_t e4x4seed;
      Float_t e2x5maxseed;
      Float_t e2x5topseed;
      Float_t e2x5bottomseed;
      Float_t e2x5leftseed;
      Float_t e2x5rightseed;      
      Float_t xseedseed;
      Float_t yseedseed;
      Float_t zseedseed;
      UInt_t  nhitsseed;
      
      //quantities from second basic cluster, if present
      Float_t ebc2;
      Float_t etabc2;
      Float_t phibc2;
      Int_t   ietabc2;
      Int_t   iphibc2;
      Int_t   ixbc2;
      Int_t   iybc2;
      Float_t etacrybc2;
      Float_t phicrybc2;
      Float_t xcrybc2;
      Float_t ycrybc2;
      Float_t thetaaxisbc2;
      Float_t phiaxisbc2;
      Float_t sigietaietabc2;
      Float_t sigiphiphibc2;
      Float_t covietaiphibc2;
      Float_t e3x3bc2;
      Float_t e5x5bc2;            
      Float_t emaxbc2;
      Float_t e2ndbc2;
      Float_t etopbc2;
      Float_t ebottombc2;
      Float_t eleftbc2;
      Float_t erightbc2;
      Float_t e1x3bc2;
      Float_t e3x1bc2;
      Float_t e1x5bc2;
      Float_t e2x2bc2;
      Float_t e4x4bc2;
      Float_t e2x5maxbc2;
      Float_t e2x5topbc2;
      Float_t e2x5bottombc2;
      Float_t e2x5leftbc2;
      Float_t e2x5rightbc2;      
      Float_t xbc2bc2;
      Float_t ybc2bc2;
      Float_t zbc2bc2;      
      UInt_t  nhitsbc2;
      
      //quantities from lowest energy basic cluster if present
      Float_t ebclast;
      Float_t etabclast;
      Float_t phibclast;
      Int_t   ietabclast;
      Int_t   iphibclast;
      Int_t   ixbclast;
      Int_t   iybclast;
      Float_t etacrybclast;
      Float_t phicrybclast;
      Float_t xcrybclast;
      Float_t ycrybclast;
      Float_t thetaaxisbclast;
      Float_t phiaxisbclast;
      Float_t sigietaietabclast;
      Float_t sigiphiphibclast;
      Float_t covietaiphibclast;
      Float_t e3x3bclast;
      Float_t e5x5bclast;             
      UInt_t  nhitsbclast;
      
      //quantities from second lowest energy basic cluster if present
      Float_t ebclast2;
      Float_t etabclast2;
      Float_t phibclast2;
      Int_t   ietabclast2;
      Int_t   iphibclast2;
      Int_t   ixbclast2;
      Int_t   iybclast2;
      Float_t etacrybclast2;
      Float_t phicrybclast2;
      Float_t xcrybclast2;
      Float_t ycrybclast2;
      Float_t thetaaxisbclast2;
      Float_t phiaxisbclast2;
      Float_t sigietaietabclast2;
      Float_t sigiphiphibclast2;
      Float_t covietaiphibclast2;
      Float_t e3x3bclast2;
      Float_t e5x5bclast2;
      UInt_t  nhitsbclast2;             
      
      //energy correction quantities from PhotonFix
      Float_t ecor;
      Float_t ecorerr;
      Float_t ecorele;
      Float_t ecoreleerr;
      Float_t etac;
      Float_t etas;
      Float_t etam;
      Float_t phic;
      Float_t phis;
      Float_t phim;  
      Float_t xz;
      Float_t xc;
      Float_t xs;
      Float_t xm;
      Float_t yz;
      Float_t yc;
      Float_t ys;
      Float_t ym;        
      
      //conversion quantities
      UChar_t hasconversion;
      Float_t convp;
      Float_t convpt;
      Float_t conveta;
      Float_t convphi;
      Float_t convdeta;
      Float_t convdphi;
      Float_t convvtxrho;
      Float_t convvtxz;
      Float_t convvtxphi;
      Float_t convleadpt;
      Float_t convtrailpt;
      Float_t convleadtrackpt;
      Char_t  convleadtrackalgo;
      Char_t  convleadtrackalgos;
      Char_t  convleadtrackcharge;
      Float_t convtrailtrackpt;
      Char_t  convtrailtrackalgo;
      Char_t  convtrailtrackalgos;      
      Char_t  trailtrackcharge;
      
      //electron quantities
      UChar_t haselectron;
      UChar_t eleisecaldriven;      
      UChar_t eleistrackerdriven;     
      Float_t elee;
      Float_t elept;
      Float_t eleeta;
      Float_t elephi;
      Char_t  elecharge;
      Float_t elefbrem;
      Float_t eledeta;
      Float_t eledphi;
      Float_t elep;
      Float_t elepin;
      Float_t elepout;
      
      //pf supercluster quantities
      UChar_t haspfsc;
      Float_t pfsce;
      Float_t pfscrawe;
      Float_t pfsceta;
      Float_t pfscphi;
      UInt_t  pfscnclusters;
      UInt_t  pfscnhits;
      Float_t pfscetawidth;
      Float_t pfscphiwidth;
      UInt_t  pfscnpsclusters;
      
      //generator level quantities
      UChar_t ispromptgen;
      Float_t gene;
      Float_t genpt;
      Float_t geneta;
      Float_t genphi;
      Float_t genz;
      Int_t   pdgid;
      Int_t   motherpdgid;
      
    // -----------------------------------------------------
    // PF-CiC4 Debug Stuff
    Float_t pfcic4_tIso1;
    Float_t pfcic4_tIso2;
    Float_t pfcic4_tIso3;
    Float_t pfcic4_covIEtaIEta;
    Float_t pfcic4_HoE;
    Float_t pfcic4_R9;
    Float_t pfcic4_wVtxInd;
    Float_t pfcic4_ecalIso3;
    Float_t pfcic4_ecalIso4;
    Float_t pfcic4_trackIsoSel03;
    Float_t pfcic4_trackIsoWorst04;
    Float_t pfcic4_combIso1;
    Float_t pfcic4_combIso2;
    // -----------------------------------------------------

    // -----------------------------------------------------
    //id mva
    //2011
    Float_t idmva_tIso1abs;
    Float_t idmva_tIso2abs;
    Float_t idmva_tIso3abs;
    Float_t idmva_absIsoEcal;
    Float_t idmva_absIsoHcal;
    //2012
    Float_t idmva_CoviEtaiPhi;
    Float_t idmva_s4ratio;
    Float_t idmva_GammaIso;
    Float_t idmva_ChargedIso_selvtx;
    Float_t idmva_ChargedIso_0p2_selvtx;
    Float_t idmva_ChargedIso_worstvtx;
    Float_t idmva_PsEffWidthSigmaRR;
    
        Float_t ebcs[NClus];
      Float_t etabcs[NClus];
      Float_t phibcs[NClus];
      Int_t ietabcs[NClus];
      Int_t iphibcs[NClus];
      Int_t ixbcs[NClus];
      Int_t iybcs[NClus];
      Float_t etacrybcs[NClus];
      Float_t phicrybcs[NClus];
      Float_t xcrybcs[NClus];
      Float_t ycrybcs[NClus];
      Float_t sigietaietabcs[NClus];
      Float_t sigiphiphibcs[NClus];
      Float_t covietaiphibcs[NClus];
      Float_t sigetaetabcs[NClus];
      Float_t sigphiphibcs[NClus];
      Float_t covetaphibcs[NClus];      
      Float_t e3x3bcs[NClus];
      Float_t e5x5bcs[NClus];            
      Float_t emaxbcs[NClus];
      Float_t e2ndbcs[NClus];
      Float_t etopbcs[NClus];
      Float_t ebottombcs[NClus];
      Float_t eleftbcs[NClus];
      Float_t erightbcs[NClus];
      Float_t e1x3bcs[NClus];
      Float_t e3x1bcs[NClus];
      Float_t e1x5bcs[NClus];
      Float_t e2x2bcs[NClus];
      Float_t e4x4bcs[NClus];
      Float_t e2x5maxbcs[NClus];
      Float_t e2x5topbcs[NClus];
      Float_t e2x5bottombcs[NClus];
      Float_t e2x5leftbcs[NClus];
      Float_t e2x5rightbcs[NClus];          
      UInt_t  nhitsbcs[NClus];      
      
      Float_t epfbcs[NClus];
      Float_t etapfbcs[NClus];
      Float_t phipfbcs[NClus];
      Int_t ietapfbcs[NClus];
      Int_t iphipfbcs[NClus];
      Int_t ixpfbcs[NClus];
      Int_t iypfbcs[NClus];
      Float_t etacrypfbcs[NClus];
      Float_t phicrypfbcs[NClus];
      Float_t xcrypfbcs[NClus];
      Float_t ycrypfbcs[NClus];
      Float_t sigietaietapfbcs[NClus];
      Float_t sigiphiphipfbcs[NClus];
      Float_t covietaiphipfbcs[NClus];
      Float_t sigetaetapfbcs[NClus];
      Float_t sigphiphipfbcs[NClus];
      Float_t covetaphipfbcs[NClus];      
      Float_t e3x3pfbcs[NClus];
      Float_t e5x5pfbcs[NClus];            
      Float_t emaxpfbcs[NClus];
      Float_t e2ndpfbcs[NClus];
      Float_t etoppfbcs[NClus];
      Float_t ebottompfbcs[NClus];
      Float_t eleftpfbcs[NClus];
      Float_t erightpfbcs[NClus];
      Float_t e1x3pfbcs[NClus];
      Float_t e3x1pfbcs[NClus];
      Float_t e1x5pfbcs[NClus];
      Float_t e2x2pfbcs[NClus];
      Float_t e4x4pfbcs[NClus];
      Float_t e2x5maxpfbcs[NClus];
      Float_t e2x5toppfbcs[NClus];
      Float_t e2x5bottompfbcs[NClus];
      Float_t e2x5leftpfbcs[NClus];
      Float_t e2x5rightpfbcs[NClus];          
      UInt_t  nhitspfbcs[NClus];      
      
      Float_t epsc[100];
      Float_t etapsc[100];
      Float_t phipsc[100];
      UChar_t planepsc[100];
      
      
  };
  
  class PhotonTreeWriterDiphotonEvent
  {
    public:
      // ------------ BTAG STUFF -------------------
      Float_t btagJet1;
      Float_t btagJet1Pt;
      Float_t btagJet1Eta;
      Float_t btagJet2;
      Float_t btagJet2Pt;
      Float_t btagJet2Eta;
      // ----------- VH/LEPTON TAG STUFF -------------
      Int_t leptonTag; // flavor-based, Moriond 2013
      Int_t VHLepTag; // MET-based, since legacy paper 2013
      Int_t VHHadTag; // VH hadronic tag
      // ----------- VERTEX SYNCHING STUFF -------

      Int_t vtxInd1;
      Int_t vtxInd2;
      Int_t vtxInd3;
	
      Float_t vtxBestPtbal  ;
      Float_t vtxBestPtasym ;
      Float_t vtxBestSumpt2 ;
      Float_t vtxBestP2Conv ;
      
      Float_t vtxMva1Z      ;
      Float_t vtxMva2Z      ;
      Float_t vtxMva3Z      ;

      Float_t vtxMva1 ;
      Float_t vtxMva2 ;
      Float_t vtxMva3 ;
      
      Int_t vtxNleg1 ;
      Int_t vtxNleg2 ;
      Int_t vtxConvIdx1 ;
      Int_t vtxConvIdx2 ;
      Int_t vtxNconv ;

      Float_t vtxConv1Z    ;
      Float_t vtxConv1DZ   ;
      Float_t vtxConv1Prob ;

      Float_t vtxConv2Z    ;
      Float_t vtxConv2DZ   ;
      Float_t vtxConv2Prob ;


      // ---------- MUON STUFF --------------------
      Float_t muonPt;
      Float_t muonEta;
      Float_t muonPhi;
      Float_t muDR1;
      Float_t muDR2;
      Float_t muIso1;
      Float_t muIso2;
      Float_t muIso3;
      Float_t muIso4;
      Float_t muD0;
      Float_t muDZ;
      Int_t muNhits;
      Float_t muChi2;
      Int_t muNpixhits;
      Int_t muNegs;
      Int_t muNMatch;
      // ------------ DIMUON STUFF ---------------
      Float_t mu1Pt;
      Float_t mu1Eta;
      Float_t mu1Phi;
      Float_t mu2Pt;
      Float_t mu2Eta;
      Float_t mu2Phi;

      // ----------- ELECTRON STUFF --------------
      Float_t elePt;
      Float_t eleEta;
      Float_t elePhi;
      Float_t eleSCEta;     
      Float_t eleIso1;
      Float_t eleIso2;
      Float_t eleIso3;
      Float_t eleIso4;
      Float_t eleDist;
      Float_t eleDcot;
      Float_t eleCoviee;
      Float_t eleDphiin;
      Float_t eleDetain;
      Float_t eleDR1;
      Float_t eleDR2;
      Float_t eleMass1;
      Float_t eleMass2;
      Int_t eleNinnerHits;     
      Float_t eleIdMva;
      // ------------ DIELECTRON STUFF -------------
      Float_t ele1Pt;
      Float_t ele1Eta;
      Float_t ele1Phi;
      Float_t ele2Pt;
      Float_t ele2Eta;
      Float_t ele2Phi;
      // -----------------------------------------
      Float_t rho;
      Float_t rho25;
      Float_t rhoold;
      Float_t genHiggspt;
      Float_t genHiggsZ;
      Float_t genmass;
      Float_t gencostheta;
      Float_t genz;
      Float_t bsX;
      Float_t bsY;
      Float_t bsZ;
      Float_t bsSigmaZ; 
      Float_t vtxX;
      Float_t vtxY;      
      Float_t vtxZ;
      Int_t   nVtx;
      Int_t   numPU;
      Int_t   numPUminus;
      Int_t   numPUplus;
      Float_t mass;
      Float_t masserr;
      Float_t masserrsmeared;
      Float_t masserrwrongvtx;
      Float_t masserrsmearedwrongvtx;
      Float_t vtxprob;
      Float_t deltamvtx;
      Float_t ptgg;
      Float_t etagg;
      Float_t phigg;
      Float_t costheta;
      Float_t massele;
      Float_t ptee;
      Float_t costhetaele;
      Float_t mt;
      Float_t cosphimet;
      Float_t mtele;
      Float_t cosphimetele;
      UInt_t  evt;
      UInt_t  run;
      UInt_t  lumi;
      UChar_t evtcat;
      UInt_t  nobj;
      Double_t fromZ; //added: Heng 2/14/2012
      Double_t fromW; //added: Heng
      Float_t zpt; //added: Heng
      Float_t allZpt; //added: Heng
      Float_t zEta; //added: Heng
      Float_t allZEta; //added: Heng

      //Met stuff
      Float_t corrpfmet;
      Float_t corrpfmetphi;
      Float_t corrpfmetx;
      Float_t corrpfmety;
      Double_t dphiMetgg;
      Double_t cosdphiMetgg;
      Double_t dphiPhPh;
      Float_t pfmet;
      Float_t pfmetphi;
      Float_t pfmetx;
      Float_t pfmety;

      
      Float_t mvametsel;
      Float_t mvametselphi;
      Float_t mvametselx;
      Float_t mvametsely;
      Float_t mvametselsig;
      
      Float_t mvametfirst;
      Float_t mvametfirstphi;
      Float_t mvametfirstx;
      Float_t mvametfirsty;      
      Float_t mvametfirstsig;
      
      Double_t spfMet;

      // ----------- VH HAD TAG STUFF -----------
      Float_t costhetastar;

      // ----------- VBF TAG STUFF -------------
      Int_t vbfTag;
      Float_t vbfbdt;

      // ----------- TTH TAG STUFF -------------
      Int_t tthTag;
      Int_t numJets;
      Int_t numBJets;
      Float_t bjetcsv;
      
      // ----------------------------------------
      UChar_t  ismc;
      Int_t mcprocid;
      Double_t mcweight;
      
      //corrected quantities from PhotonFix corrections
      Float_t masscor;
      Float_t masscorerr;
      Float_t masscorele;
      Float_t masscoreleerr;

      //jet quantities for met phi cut
      Float_t jetleadNoIDpt;
      Float_t jetleadNoIDeta;
      Float_t jetleadNoIDphi;
      Float_t jetleadNoIDmass;

      //jet quantities
      Float_t jet1pt;
      Float_t jet1eta;
      Float_t jet1phi;
      Float_t jet1mass;
      Float_t jet2pt;
      Float_t jet2eta;
      Float_t jet2phi;
      Float_t jet2mass;
      Float_t jetcentralpt;
      Float_t jetcentraleta;
      Float_t jetcentralphi;
      Float_t jetcentralmass;
      Float_t dijetpt;
      Float_t dijeteta;
      Float_t dijetphi;
      Float_t dijetmass;
      Float_t jetetaplus;
      Float_t jetetaminus;
      
      //      Float_t uncorrjet1pt;
      //      Float_t uncorrjet1eta;
      //      Float_t uncorrjet1phi;
      //      Float_t uncorrjet1mass;
      //      Float_t uncorrjet2pt;
      //      Float_t uncorrjet2eta;
      //      Float_t uncorrjet2phi;
      //      Float_t uncorrjet2mass;
      //      Float_t uncorrjetcentralpt;
      //      Float_t uncorrjetcentraleta;
      //      Float_t uncorrjetcentralphi;
      //      Float_t uncorrjetcentralmass;
      //      Float_t diuncorrjetpt;
      //      Float_t diuncorrjeteta;
      //      Float_t diuncorrjetphi;
      //      Float_t diuncorrjetmass;
      //      Float_t uncorrjetetaplus;
      //      Float_t uncorrjetetaminus;

      Float_t zeppenfeld;
      Float_t dphidijetgg;
      
      PhotonTreeWriterPhoton<16> photons[2];
  }; // PhotonTreeWriterDiphotonEvent
  
  
  class PhotonTreeWriter : public BaseMod
  {
  public:
    enum CounterFlag {kCountNotDefined = -1};
    typedef std::pair<const Particle*, double>  DeltaRVeto;
    typedef std::vector<DeltaRVeto>             DeltaRVetoVector;
    typedef std::vector<const PFJet*>           PFJetVector;
    PhotonTreeWriter(const char *name ="PhotonTreeWriter",
		     const char *title="Selecting PhotonPairs");
    
    ~PhotonTreeWriter();

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
    void                SetPFSuperClusterName(const char *n) { fPFSuperClusterName = n;      }
    void                SetPFJetName(const char *n)       { fPFJetName = n;              }
    void                SetGenJetName(const char *n)      { fGenJetName = n;             }
    void                SetuncorrPFJetName(const char *n) { funcorrPFJetName = n;        }
    void                SetPFNoPileUpName(const char *n)  { fPFNoPileUpName  = n;        }
    void                SetPFPileUpName(const char *n)    { fPFPileUpName  = n;          }


    void                SetPFJetsFromBranch(Bool_t b)     { fPFJetsFromBranch = b;       }
    void                SetEnableJets(Bool_t b)           { fEnableJets = b;             }
    void                SetEnableGenJets(Bool_t b)        { fEnableGenJets = b;          }
    void                SetApplyJetId(Bool_t b)           { fApplyJetId = b;             }
    void                SetApplyLeptonTag(Bool_t b)       { fApplyLeptonTag = b;         }
    void                SetApplyVHLepTag(Bool_t b)        { fApplyVHLepTag = b;          }
    void                SetApplyVHHadTag(Bool_t b)        { fApplyVHHadTag = b;          }
    void                SetApplyVBFTag(Bool_t b)          { fApplyVBFTag = b;            }
    void                SetApplyTTHTag(Bool_t b)          { fApplyTTHTag = b;            }
    void                SetApplyBTag(Bool_t b)            { fApplyBTag = b;              }
    void                SetApplyPFMetCorr(Bool_t b)       { fApplyPFMetCorrections = b;  }
    void                SetPhFixDataFile(const char *n)   { fPhFixDataFile = n;          }
    void                SetVerbosityLevel(Int_t b)       { fVerbosityLevel = b;         }




    // set basic Cut variables (FOR PRE-SELECTION)

    // is Data Or Not?
    void                SetIsData     (Bool_t b)             { fIsData = b; };
    void                SetIsCutBased (Bool_t b)             { fIsCutBased = b; };
    

    void                SetApplyElectronVeto(Bool_t b)   { fApplyElectronVeto = b;     }          

    void                SetTupleName(const char* c)          { fTupleName = c; }
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name)    { fGoodElectronName = name; }
    void                SetWriteDiphotonTree(Bool_t b)       { fWriteDiphotonTree = b; }
    void                SetWriteSingleTree(Bool_t b)         { fWriteSingleTree = b; }
    void                SetLoopOnGoodElectrons(Bool_t b)     { fLoopOnGoodElectrons = b; }
    void                SetEnablePFPhotons(Bool_t b)         { fEnablePFPhotons = b; }
    void                SetExcludeSinglePrompt(Bool_t b)     { fExcludeSinglePrompt = b; }
    void                SetExcludeDoublePrompt(Bool_t b)     { fExcludeDoublePrompt = b; }

    void                SetLeptonTagElectronsName(TString name) { fLeptonTagElectronsName = name; }
    void                SetLeptonTagMuonsName    (TString name) { fLeptonTagMuonsName     = name; }
    void                SetLeptonTagSoftElectronsName(TString name) { fLeptonTagSoftElectronsName = name; }
    void                SetLeptonTagSoftMuonsName    (TString name) { fLeptonTagSoftMuonsName     = name; }
    void                SetFillClusterArrays(Bool_t b)          { fFillClusterArrays = b; }

    void                SetDo2012LepTag(Bool_t b) {                         fDo2012LepTag = b; }
    
    void                SetBeamspotWidth(Double_t x)            { fBeamspotWidth = x; }

    void                SetElectronMVAWeightsSubdet0Pt10To20(TString s)
                        { fElectronMVAWeights_Subdet0Pt10To20  = s; }
    void                SetElectronMVAWeightsSubdet1Pt10To20(TString s)
                        { fElectronMVAWeights_Subdet1Pt10To20  = s; }
    void                SetElectronMVAWeightsSubdet2Pt10To20(TString s)
                        { fElectronMVAWeights_Subdet2Pt10To20  = s; }
    void                SetElectronMVAWeightsSubdet0Pt20ToInf(TString s)
                        { fElectronMVAWeights_Subdet0Pt20ToInf = s; }
    void                SetElectronMVAWeightsSubdet1Pt20ToInf(TString s)
                        { fElectronMVAWeights_Subdet1Pt20ToInf = s; }
    void                SetElectronMVAWeightsSubdet2Pt20ToInf(TString s)
                        { fElectronMVAWeights_Subdet2Pt20ToInf = s; }

    void                SetDoSynching(bool b) {fDoSynching = b;}

    template<typename Element, typename Collection>
    static UInt_t       IndexOfElementInCollection(
                          const Element *element,
                          const Collection *collection
                          );

    static UInt_t       IndexOfNearestSuperClusterInCollection(
                          const SuperCluster    *element,
                          const SuperClusterCol *collection
                          );

    static double       SuperClusterPairMass(
                          const SuperCluster *sc1,
                          const SuperCluster *sc2
                          );

    static FourVectorM  SuperClusterFourVectorM(const SuperCluster *sc1);

    void                SetR9Rescale(bool dor9rescale, double p0b, double p1b, double  p0e,  double  p1e) {fdor9rescale = dor9rescale; fp0b = p0b; fp1b = p1b; fp0e = p0e; fp1e = p1e;}

  protected:
    void                Process();
    void                SlaveBegin();
    void                SlaveTerminate();
    void                Terminate();
    // Private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z, Float_t& mass);
    Float_t             GetEventCat    (PhotonTools::CiCBaseLineCats cat1,
					PhotonTools::CiCBaseLineCats cat2);
    void                ApplyLeptonTag(const Photon *phHard,
                                       const Photon *phSoft,
                                       const Vertex *selvtx);
    const Muon *        GetLeptonTagMuon(const Photon *phHard,
                                         const Photon *phSoft);
    void                SetLeptonTagMuonVars(const Photon *phHard,
                                             const Photon *phSoft,
                                             const Muon *muon);
    const Electron *    GetLeptonTagElectron(const Photon *phHard,
                                             const Photon *phSoft);
    bool                MassOfPairIsWithinWindowAroundMZ(
                          const Particle * particle1,
                          const Particle * particle2,
                          Float_t halfWindowSize,
                          Float_t MZ = 91.19
                        );
    void                ApplyVHLepTag(const Photon *phHard,
                                      const Photon *phSoft,
                                      const Vertex *selvtx);
    void                ApplyVHHadTag(const Photon *phHard,
                                      const Photon *phSoft,
                                      const Vertex *selvtx,
                                      const Jet    *jet1,
                                      const Jet    *jet2);
    void                ApplyTTHTag(const Photon *phHard, 
                                    const Photon *phSoft,
                                    const Vertex *selvtx);
    bool                VHLepHasDielectron(const Photon *phHard,
                                           const Photon *phSoft);
    bool                VHLepHasDimuon(const Photon *phHard,
                                      const Photon *phSoft);
    UInt_t              VHLepNumberOfJets(const Photon *phHard,
                                          const Photon *phSoft,
                                          const Vertex *selvtx,
                                          const Particle *lepton);
    UInt_t              NumberOfJets(const Photon *phHard,
                                     const Photon *phSoft,
                                     const Vertex *selvtx,
                                     const double minJetPt,
                                     const double maxAbsEta);
    UInt_t              NumberOfBJets(const Photon *phHard,
                                      const Photon *phSoft,
                                      const Vertex *selvtx,
                                      const double minJetPt,
                                      const double maxAbsEta);
    PFJetVector *      GetSelectedPFJets(const DeltaRVetoVector &drVetos,
                                         const Vertex &vertex,
                                         const double minJetPt,
                                         const double maxJetAbsEta);
    PFJetVector *      GetSelectedPFBJets(const PFJetVector &pfjets);
    const Particle *    TTHSelectLepton(const Photon *phHard,
                                        const Photon *phSoft,
                                        const Vertex *selvtx);
    const Electron *    TTHSelectElectron(const Photon *phHard,
                                          const Photon *phSoft,
                                          const Vertex *selvtx);
    const Muon *        TTHSelectMuon(const Photon *phHard,
                                      const Photon *phSoft,
                                      const Vertex *selvtx);
    void                LogEventInfo();
    double              GetElectronIdMva(const Electron *electron);
    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fPFPhotonName;
    TString             fElectronName;
    TString             fGoodElectronName;
    TString             fConversionName;
    TString             fPFConversionName;
    TString             fTrackBranchName;
    TString             fPileUpDenName;
    TString             fPVName;
    TString             fBeamspotName;
    TString             fPFCandName;
    TString             fPFNoPileUpName;  //name of pfnpu collection
    TString             fPFPileUpName;    //name of pfpu collection

    TString             fMCParticleName;
    TString             fMCEventInfoName;

    TString             fPileUpName;
    TString             fPFSuperClusterName;
    TString             fPFMetName;
    TString             fPFJetName;

    TString             funcorrPFJetName;
    TString             fGenJetName;   //added to do pfmet correction 05/01/2012

    TString             fLeptonTagElectronsName;
    TString             fLeptonTagMuonsName;
    TString             fLeptonTagSoftElectronsName;
    TString             fLeptonTagSoftMuonsName;

    
    // is it Data or MC?
    Bool_t              fIsData;
    // Is it cut-based? "false" means mass-factorized.
    Bool_t              fIsCutBased;
    
    // in case there's some PV pre-selection
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;
    Bool_t              fPFJetsFromBranch;

    Bool_t              fDoSynching;

    const PhotonCol               *fPhotons;
    const PhotonCol               *fPFPhotons;
    const ElectronCol             *fElectrons;
    const ElectronCol             *fGoodElectrons;    
    const DecayParticleCol        *fConversions;
    const DecayParticleCol        *fPFConversions;
    const TrackCol                *fTracks;
    const PileupEnergyDensityCol  *fPileUpDen;
    const VertexCol               *fPV;
    const BeamSpotCol             *fBeamspot;
    const PFCandidateCol          *fPFCands;
    const MCParticleCol           *fMCParticles;
    const MCEventInfo             *fMCEventInfo;
    const PileupInfoCol           *fPileUp;    
    const SuperClusterCol         *fPFSuperClusters;   
    const PFMetCol                *fPFMet;
    const JetCol                  *fPFJets;
    const GenJetCol               *fGenJets;
    const PFJetCol                *funcorrPFJets;
    
    const ElectronCol             *fLeptonTagElectrons;
    const MuonCol                 *fLeptonTagMuons;
    const ElectronCol             *fLeptonTagSoftElectrons;
    const MuonCol                 *fLeptonTagSoftMuons;
    const PFCandidateCol          *fPFNoPileUpCands;  //!pfnpu collection
    const PFCandidateCol          *fPFPileUpCands;    //!pfpu collection

    // --------------------------------
    Bool_t                         fLoopOnGoodElectrons; //loop over good elecs instead of photons
    Bool_t                         fApplyElectronVeto;   //invert elec veto (for cic sel. only atm)
    Bool_t                         fWriteDiphotonTree;
    Bool_t                         fWriteSingleTree;

    Bool_t                         fEnablePFPhotons;
    
    Bool_t                         fExcludeSinglePrompt;
    Bool_t                         fExcludeDoublePrompt;
    
    Bool_t                         fEnableJets;
    Bool_t                         fEnableGenJets;
    Bool_t                         fApplyJetId;

    Bool_t                         fApplyLeptonTag;
    Bool_t                         fApplyVHLepTag;
    Bool_t                         fApplyVHHadTag;
    Bool_t                         fApplyVBFTag;
    Bool_t                         fApplyTTHTag;
    Bool_t                         fApplyBTag;
    Bool_t                         fApplyPFMetCorrections;

    Bool_t                         fFillClusterArrays;
    Bool_t                         fFillVertexTree;
    
    Bool_t                         fDo2012LepTag;
    
    Int_t                          fVerbosityLevel;

    TString                        fPhFixDataFile;
    PhotonFix                      fPhfixph;
    PhotonFix                      fPhfixele;
    
    Double_t                       fBeamspotWidth;
    
    TFile                          *fTmpFile;

    // --------------------------------
    // variables for vbf
    float jet1pt_vbf;
    float jet2pt_vbf;
    float deltajeteta_vbf;
    float dijetmass_vbf;
    float zeppenfeld_vbf;
    float dphidijetgg_vbf;
    float diphoptOverdiphomass_vbf;
    float pho1ptOverdiphomass_vbf;
    float pho2ptOverdiphomass_vbf; 

    // --------------------------------
    // validation Tuple
    TString                        fTupleName;
    PhotonTreeWriterDiphotonEvent* fDiphotonEvent;
    PhotonTreeWriterPhoton<16>*    fSinglePhoton;    
    TTree*                         hCiCTuple;
    TTree*                         hCiCTupleSingle;

    PhotonTreeWriterVtx           *fDiphotonVtx;
    TTree                         *hVtxTree;
    
    MVAMet                         fMVAMet;
    JetIDMVA                       fJetId;
    MVAVBF                         fMVAVBF;   
    
    VertexTools           fVtxTools;

    ElectronIDMVA                 *fElectronIDMVA;
    TString                   fElectronMVAWeights_Subdet0Pt10To20;
    TString                   fElectronMVAWeights_Subdet1Pt10To20;
    TString                   fElectronMVAWeights_Subdet2Pt10To20;
    TString                   fElectronMVAWeights_Subdet0Pt20ToInf;
    TString                   fElectronMVAWeights_Subdet1Pt20ToInf;
    TString                   fElectronMVAWeights_Subdet2Pt20ToInf;

    //R9 rescale
    bool                  fdor9rescale;                   
    double                fp0b;                          
    double                fp1b;                          
    double                fp0e;                          
    double                fp1e;     

    RhoUtilities::RhoType    fTheRhoType;
    UInt_t                         fProcessedEvents;

    ClassDef(PhotonTreeWriter, 1) // Photon identification module
      };
}
#endif
