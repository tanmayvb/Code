// -*- C++ -*-
//
// Package:    Test/Triggereffi
// Class:      Triggereffi
// 
/**\class Triggereffi Triggereffi.cc Test/Triggereffi/plugins/Triggereffi.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tanmay Sarkar
//         Created:  Wed, 11 Jan 2017 06:26:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <map>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TFormula.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <cmath>
#include "TMath.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TRandom.h"

#include "TH2F.h"
#include "TProfile.h"
#include <fstream>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReport.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReportEntry.h"
#include "CondFormats/DataRecord/interface/L1GtStableParametersRcd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

//#include "Test/QCDEventShape/plugins/EventShape_vector.h"


#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <CondFormats/DataRecord/interface/JetResolutionRcd.h>
#include <CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h>
#include "FWCore/Utilities/interface/typelookup.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

using namespace edm;
using namespace reco;
using namespace std;
using namespace CLHEP;
using namespace trigger;
using namespace math;

const int nHLTmx=8;
const int njetetamn=1;
const char* jethlt_name[nHLTmx]={"HLT_DiPFJetAve60_v","HLT_DiPFJetAve80_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve320_v", "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve500_v"};
double l1Pt[nHLTmx] = {0,35,60,90,120,170,170,170};
double jethlt_thr[nHLTmx]={60,80,140,200,260,320,400,500};
double leadingPtThreshold[nHLTmx]={60,80,140,200,260,320,400,500};
const char* jethlt_lowest={"HLT_DiPFJetAve40_v"};
double etarange[njetetamn] ={2.4};
bool trgpas[nHLTmx];//={0,0,0,0,0,0,0,0};
// class declaration
//
int l1pres[nHLTmx], hltpres[nHLTmx], compres[nHLTmx];
struct triggervar{
  HepLorentzVector trg4v;
  bool            both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
};
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
//

class Triggereffi : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Triggereffi(const edm::ParameterSet&);
      ~Triggereffi();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //    TH1F* hlt_dijet[nHLTmx];
  //    TH1F* hlt_dijetref[nHLTmx];
  
  //Dijet trigger efficiency
  TH1F* hlt_dijettag[nHLTmx][njetetamn];
  TH1F* hlt_dijetprob[nHLTmx][njetetamn];
  
  std::string theRootFileName;
  std::string theHLTTag;
  // ----------member data ---------------------------
  edm::EDGetTokenT<GenEventInfoProduct> generator1_;
  edm::EDGetTokenT<pat::JetCollection> jetSrcToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genSrcToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> PFSrcToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::EDGetTokenT<reco::GenJetCollection> genjetToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileup_;
  edm::EDGetTokenT<reco::PFJetCollection> ak5PFjetToken_;
  edm::EDGetTokenT<reco::GenJetCollection> ak5GenJetToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfCTEQWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfMMTHWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfNNPDFWeightsInputToken_;
  const edm::EDGetTokenT<LHERunInfoProduct> LHERunInfoToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
  edm::EDGetTokenT<double> m_rho_token;
   HLTPrescaleProvider hltPrescaleProvider_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Triggereffi::Triggereffi(const edm::ParameterSet& iConfig):
  generator1_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("evtinfo"))),
  jetSrcToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  genSrcToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genSrc"))),
  PFSrcToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bsSrc"))),
  genjetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genjetSrc"))),
  pileup_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSrc"))),
  ak5PFjetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak5pfJetSrc"))),
  ak5GenJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak5genJetSrc"))),
  pdfCTEQWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFCTEQWeightsInputTag"))),
  pdfMMTHWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFMMTHWeightsInputTag"))),
  pdfNNPDFWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFNNPDFWeightsInputTag"))),
  LHERunInfoToken_(consumes<LHERunInfoProduct, edm::InRun >(iConfig.getParameter<edm::InputTag>("LHERunInfoProductInputTag"))),
  lheEventProductToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductInputTag"))),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)

{
   //now do what ever initialization is needed
   //usesResource("TFileService");
  // theRootFileName = iConfig.getUntrackedParameter<string>("RootFileName");
   edm::Service<TFileService> fs;
   theHLTTag = iConfig.getUntrackedParameter<string>("HLTTag", "HLT");
   char name[200];
   char title[200];
   
   
   /*   for (int ij=0; ij<nHLTmx; ij++) {
   // for (int jk=0; jk<1; jk++) {
   sprintf(name, "hlt_dijetref_%i", ij);
   sprintf(title, "dijetref probed P_T : (%s) |i#eta|<", jethlt_name[ij]);
   hlt_dijetref[ij] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
   hlt_dijetref[ij]->Sumw2();
   sprintf(name, "hlt_dijetprob_%i", ij);
   sprintf(title, "dijet probed P_T : (%s) |i#eta|<", jethlt_name[ij]);
   hlt_dijet[ij] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
   hlt_dijet[ij]->Sumw2();
   //}
   }*/
   for (int ij=0; ij<nHLTmx; ij++) {
     for (int jk=0; jk<njetetamn; jk++) {
       sprintf(name, "hlt_dijettag_%i_%i", ij, jk);
       sprintf(title, "dijet tagged P_T : (%s) |i#eta|<%g", jethlt_name[ij], etarange[jk]);
       hlt_dijettag[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
       hlt_dijettag[ij][jk]->Sumw2();
       
       sprintf(name, "hlt_dijetprob_%i_%i", ij, jk);
       sprintf(title, "dijet probed P_T : (%s) |i#eta|<%g", jethlt_name[ij], etarange[jk]);
       hlt_dijetprob[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
       hlt_dijetprob[ij][jk]->Sumw2();
     }
   }
   
   
   
}


Triggereffi::~Triggereffi()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Triggereffi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/
const char* variab1;
//const char* variab2;
 double aveleadingpt =0;
 edm::Handle<edm::TriggerResults> trigRes;
 iEvent.getByToken(triggerBits_, trigRes);
 
 edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
 iEvent.getByToken(triggerObjects_, triggerObjects);
 
 edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
 iEvent.getByToken(triggerPrescales_, triggerPrescales);
 const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
 bool isInEtaRange[njetetamn]={0};
 edm::Handle<pat::JetCollection> ak4PFJets;
 iEvent.getByToken(jetSrcToken_, ak4PFJets);
 vector<triggervar> alltrgobj;
 alltrgobj.clear();
 if ((!ak4PFJets.isValid()) ||  ak4PFJets->size() <2) return;
 if (ak4PFJets.isValid() &&  ak4PFJets->size()>=2) {
    for (int iet=0; iet<njetetamn; iet++) {
      isInEtaRange[iet] = true;
    }    
    
   for (int ij=0; ij<2; ij++) {
     for (int iet=0; iet<njetetamn; iet++) {
        if (abs((*ak4PFJets)[ij].eta())>etarange[iet]) { isInEtaRange[iet] = false;}
      }
     double NHF = (*ak4PFJets)[ij].neutralHadronEnergyFraction();
     double NEMF = (*ak4PFJets)[ij].neutralEmEnergyFraction();
     double CHF = (*ak4PFJets)[ij].chargedHadronEnergyFraction();
     double CEMF = (*ak4PFJets)[ij].chargedEmEnergyFraction();
     int NumConst = (*ak4PFJets)[ij].chargedMultiplicity()+(*ak4PFJets)[ij].neutralMultiplicity();
     int NumNeutralParticles =(*ak4PFJets)[ij].neutralMultiplicity();
     int CHM = (*ak4PFJets)[ij].chargedMultiplicity();
     bool looseJetID =false;
     if(abs((*ak4PFJets)[ij].eta())<=3.0){
       if( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs((*ak4PFJets)[ij].eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) ||
						     
						     abs((*ak4PFJets)[ij].eta())>2.4) ) looseJetID =true;
     } else {
       if( (NEMF<0.90 && NumNeutralParticles>10) ) looseJetID =true;
     }
     
     if (abs((*ak4PFJets)[ij].eta())>3.0) {looseJetID = false;}
     if ((*ak4PFJets)[ij].pt()<30.0) {looseJetID = false;}
     
     if (looseJetID) {
       aveleadingpt +=(*ak4PFJets)[ij].pt();
     } else {
       aveleadingpt -=100000;
     }
   }
   aveleadingpt /=2.0;
 }
// trgpas[nHLTmx]={0,0,0,0,0,0,0,0};
 bool trg_prev=false;
 // trgpas[nHLTmx]={false};
 for (int jk=-1; jk<nHLTmx; jk++) {
   for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
     std::string name = names.triggerName(ij);
     variab1 = name.c_str();
     if ((jk<0 && strstr(variab1,jethlt_lowest) && strlen(variab1)-strlen(jethlt_lowest)<5) ||
	 (jk>=0 && strstr(variab1,jethlt_name[jk]) && strlen(variab1)-strlen(jethlt_name[jk])<5)) {
       
    //   const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltPrescaleProvider_.prescaleValuesInDetail(iEvent,iSetup,variab1));
       if (jk>=0) {
//	 l1pres[jk] = prescalesInDetail.first[0].second;
	 
//	 if (jk>=3 && l1pres[jk]>1) { l1pres[jk]=1.0;}
	 
	 
//	 hltpres[jk] = prescalesInDetail.second;
	 compres[jk] = triggerPrescales->getPrescaleForIndex(ij);
	 if (trigRes->accept(ij)) {trgpas[jk] = true;}
   //     cout<< "OK1"<< variab1 <<endl;
	 if (trg_prev){
   //     cout<< "OK2"<<endl;
	   for (int iet=0; iet<njetetamn; iet++) {
	     if (isInEtaRange[iet] && (aveleadingpt>30)) {
	       hlt_dijettag[jk][iet]->Fill(aveleadingpt,compres[jk]);
               cout << ij << " ; "<< variab1 << ";  tagpt = " << aveleadingpt << "comp scale = "<< compres[jk] << endl;
	       if (trigRes->accept(ij)) {hlt_dijetprob[jk][iet]->Fill(aveleadingpt, compres[jk]);//} //{, (isMC) ? 1.0 : compres[jk]);}
                cout << "==============" << trigRes->accept(ij) <<endl;
               cout << ij << " ; "<< variab1 << ";  probpt = " << aveleadingpt <<"comp scale = "<< compres[jk] << endl;
                cout << "==============" << endl;
              }
	     }
	   }
	 }
    //     cout<< "OK3"<<endl;
	 trg_prev = trigRes->accept(ij);
	 break;
       } else {
	 trg_prev = trigRes->accept(ij);
	 break;
       }
     }
   }
 }
 
 
 
}


// ------------ method called once each job just before starting event loop  ------------

void 
Triggereffi::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Triggereffi::endJob() 
{
}
/*
void
Triggereffi::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
cout << "Write test 4 = ok " << endl;
        bool changed(true);
  if (hltPrescaleProvider_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
  HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();
    hltConfig.dump("Triggers");
    hltConfig.dump("PrescaleTable");

       } else {
         }

}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Triggereffi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Triggereffi);
