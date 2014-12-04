/// \class l1t::Stage1Layer2TauAlgorithm
///
/// Description: interface for MP firmware
///
/// Implementation:
///
/// \author: Kalanand Mishra - Fermilab
///
/// Tau definition: 4x8 towers.


#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2TauAlgorithmImp.h"
#include "DataFormats/L1TCalorimeter/interface/CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"
#include "L1Trigger/L1TCalorimeter/interface/PUSubtractionMethods.h"
#include "L1Trigger/L1TCalorimeter/interface/JetFinderMethods.h"
#include "L1Trigger/L1TCalorimeter/interface/legacyGtHelper.h"


using namespace std;
using namespace l1t;


Stage1Layer2TauAlgorithmImpPP::Stage1Layer2TauAlgorithmImpPP(CaloParamsStage1* params) : params_(params)
{
}

Stage1Layer2TauAlgorithmImpPP::~Stage1Layer2TauAlgorithmImpPP(){};




void l1t::Stage1Layer2TauAlgorithmImpPP::processEvent(const std::vector<l1t::CaloEmCand> & EMCands,
						      const std::vector<l1t::CaloRegion> & regions,
						      const std::vector<l1t::Jet> * jets,
						      std::vector<l1t::Tau> * taus) {

  double towerLsb = params_->towerLsbSum();

  std::string regionPUSType = params_->regionPUSType();
  std::vector<double> regionPUSParams = params_->regionPUSParams();
  int tauSeedThreshold= floor( params_->tauSeedThreshold()/towerLsb + 0.5); // convert GeV to HW units
  int tauNeighbourThreshold= floor( params_->tauNeighbourThreshold()/towerLsb + 0.5); // convert GeV to HW units
  int jetSeedThreshold= floor( params_->jetSeedThreshold()/towerLsb + 0.5); // convert GeV to HW units
  int switchOffTauVeto = floor( params_->switchOffTauVeto()/towerLsb + 0.5);
  int switchOffTauIso = floor( params_->switchOffTauIso()/towerLsb + 0.5);
  double tauRelativeJetIsolationLimit = params_->tauRelativeJetIsolationLimit();
  double tauRelativeJetIsolationCut = params_->tauRelativeJetIsolationCut();

  std::vector<l1t::CaloRegion> *subRegions = new std::vector<l1t::CaloRegion>();



  //Region Correction will return uncorrected subregions if
  //regionPUSType is set to None in the config
  RegionCorrection(regions, subRegions, regionPUSParams, regionPUSType);



  // ----- need to cluster jets in order to compute jet isolation ----
  std::vector<l1t::Jet> *unCorrJets = new std::vector<l1t::Jet>();
  //slidingWindowJetFinder(jetSeedThreshold, subRegions, unCorrJets);
  TwelveByTwelveFinder(jetSeedThreshold, subRegions, unCorrJets);

  std::vector<l1t::Tau> *preGtTaus = new std::vector<l1t::Tau>();

  //std::cout << "starting regions loop" << std::endl;
  for(CaloRegionBxCollection::const_iterator region = subRegions->begin();
      region != subRegions->end(); region++) {

    int regionEt = region->hwPt();
    if(regionEt < tauSeedThreshold) continue;

    int regionEta = region->hwEta();
    int regionPhi = region->hwPhi();

    //int associatedSecondRegionEt =
    //  AssociatedSecondRegionEt(region->hwEta(), region->hwPhi(),
    //			       *subRegions);

    int tauEt=regionEt;
    int isoFlag=0;  // is 1 if it passes the relative jet iso requirement
    int quality = 1;  //doesn't really mean anything and isn't used

    int highestNeighborEt=0;
    int highestNeighborEta=999;
    int highestNeighborPhi=999;
    int highestNeighborTauVeto=999;

    // if (regionEt>0) std::cout << "CCLA Prod: TauVeto: " << region->hwQual() << "\tET: " << regionEt << "\tETA: " << regionEta  << "\tPhi: " << regionPhi  << std::endl;

    //Find neighbor with highest Et
    for(CaloRegionBxCollection::const_iterator neighbor = subRegions->begin();
	neighbor != subRegions->end(); neighbor++) {

      
      int highestSubNeighborEt=0;
      int highestSubNeighborEta=999;
      int highestSubNeighborPhi=999;
      int highestSubNeighborTauVeto=999;

      int neighborPhi = neighbor->hwPhi();
      int neighborEta = neighbor->hwEta();
      int deltaPhi = regionPhi - neighborPhi;
      if (std::abs(deltaPhi) == L1CaloRegionDetId::N_PHI-1)
	      deltaPhi = -deltaPhi/std::abs(deltaPhi); //18 regions in phi

      deltaPhi = std::abs(deltaPhi);
      int deltaEta = std::abs(regionEta - neighborEta);
      	
      if (deltaPhi + deltaEta > 0 && deltaPhi + deltaEta < 2) {  //nondiagonal neighbors
        std::cout<<"neighbor et, highest neighbor et  " << neighbor->hwPt() << " , " << highestNeighborEt<<std::endl;
	if (neighbor->hwPt() > highestNeighborEt || (highestNeighborPhi == 999)) {
          std::cout << "neighbor has highest et yet, looping through subneighbors" << std::endl;
// find neighbor of neighbor with highest Et
          for(CaloRegionBxCollection::const_iterator subneighbor = subRegions->begin();
              subneighbor != subRegions->end(); subneighbor++) {
  
            int subneighborPhi = subneighbor->hwPhi();
            int subneighborEta = subneighbor->hwEta();
            int subdeltaPhi = neighborPhi - subneighborPhi;
            if (std::abs(subdeltaPhi) == L1CaloRegionDetId::N_PHI-1)
              subdeltaPhi = -subdeltaPhi/std::abs(subdeltaPhi); //18 regions in phi
     
            subdeltaPhi = std::abs(subdeltaPhi);
            int subdeltaEta = std::abs(neighborEta - subneighborEta);
     
            if (subdeltaPhi + subdeltaEta > 0 && subdeltaPhi + subdeltaEta < 2) {  //nondiagonal subneighbors
            std::cout<<"subneighbor et, highest subneighbor et  " << subneighbor->hwPt() << " , " << highestSubNeighborEt<<std::endl;
              if (subneighbor->hwPt() > highestSubNeighborEt || highestSubNeighborPhi == 999){
                std::cout << "sub neighbor has highest et yet" << std::endl;
                highestSubNeighborEt = subneighbor->hwPt();
                highestSubNeighborEta = subneighbor->hwEta();
                highestSubNeighborPhi = subneighbor->hwPhi();
                int subneighborTauVeto = subneighbor->hwQual() & 0x1; // tauVeto should be the first bit of quality integer
                highestSubNeighborTauVeto = subneighborTauVeto;
              }
            }
          }
          // check to make sure we aren't double counting energy before declaring a highest neighbor
          string NESWNeighbor = findNESW(neighborEta, neighborPhi, highestSubNeighborEta, highestSubNeighborPhi); //NESW between neighbor and subneighbor
          string NESWSubNeighbor = findNESW(regionEta, regionPhi, highestSubNeighborEta, highestSubNeighborPhi); //NESW between central and subneighbor (overlap check)
          string NESWCandidate = findNESW(regionEta, regionPhi, neighborEta, neighborPhi); //NESW between central and neighbor
          std::cout<<"tau et, subneighbor et, NESWNeighbor, NESWSubNeighbor, NESWCandidate" << tauEt << " , " << highestSubNeighborEt << " , " << NESWNeighbor << " , " << NESWSubNeighbor  << ", " << NESWCandidate << std::endl;
          if(NESWNeighbor == "isNorth" || NESWNeighbor == "isSouth" || NESWNeighbor == "isEast" || NESWNeighbor == "isWest"){
            std::cout<<"highest subneighbor is non diagonal to highest neighbor" << std::endl;
            if(NESWSubNeighbor == "isOverlap" || (tauEt > highestSubNeighborEt && (NESWCandidate=="isEast" || NESWCandidate=="isNorth"))
              || (tauEt >= highestSubNeighborEt && (NESWCandidate=="isSouth" || NESWCandidate=="isWest"))){ //check if central is highest neighbor of neighbor before storing highest neighbor info
              std::cout << "adding highest neighbor info" << std::endl; 
   
    	      highestNeighborEt = neighbor->hwPt();
	      highestNeighborEta = neighbor->hwEta();
	      highestNeighborPhi = neighbor->hwPhi();
	      int neighborTauVeto = neighbor->hwQual() & 0x1; // tauVeto should be the first bit of quality integer
	      highestNeighborTauVeto = neighborTauVeto;
	    }
            if (highestSubNeighborEt > tauEt){
              std::cout << "should skip recording neighbor as highest neighbor" << std::endl;  
	    }
          }
        }
      }

    }
    string NESW = findNESW(regionEta, regionPhi, highestNeighborEta, highestNeighborPhi);

    std::cout << "tau et, neighbor et " << tauEt << " " << highestNeighborEt << std::endl;
    if((tauEt > highestNeighborEt && (NESW=="isEast" || NESW=="isNorth"))
       || (tauEt >= highestNeighborEt && (NESW=="isSouth" || NESW=="isWest"))
       || highestNeighborEt == 0 ) {

      if (highestNeighborEt >= tauNeighbourThreshold) tauEt += highestNeighborEt;

      int regionTauVeto = region->hwQual() & 0x1;  // tauVeto should be the first bit of quality integer
      //std::cout<< "regiontauveto, neighbor " << regionTauVeto << " " << highestNeighborTauVeto << std::endl;

      double jetIsolation = JetIsolation(tauEt, region->hwEta(), region->hwPhi(), *unCorrJets);
      if ((highestNeighborTauVeto == 0 && regionTauVeto == 0) || tauEt > switchOffTauVeto) {
	if (jetIsolation < tauRelativeJetIsolationCut || (tauEt >= switchOffTauIso && jetIsolation < tauRelativeJetIsolationLimit)
	    || (std::abs(jetIsolation - 999.) < 0.1) ) isoFlag=1;
      
      //if (highestSubNeighborEt > tauEt){
          //std::cout << "!!!!!" << std::endl;
        //std::cout << "tau et, neighbor et, subneighbor et " << tauEt << " " << highestNeighborEt << " " << highestSubNeighborEt<< std::endl;
            //}
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > tauLorentz(0,0,0,0);

        l1t::Tau theTau(*&tauLorentz, tauEt, region->hwEta(), region->hwPhi(), quality, isoFlag);
        preGtTaus->push_back(theTau);


      }
    }
  }
  
  TauToGtScales(params_, preGtTaus, taus);

  delete subRegions;
  delete unCorrJets;
  delete preGtTaus;

  //the taus should be sorted, highest pT first.
  // do not truncate the tau list, GT converter handles that
  auto comp = [&](l1t::Tau i, l1t::Tau j)-> bool {
    return (i.hwPt() < j.hwPt() );
  };

  std::sort(taus->begin(), taus->end(), comp);
  std::reverse(taus->begin(), taus->end());
}





//  Compute jet isolation.
double l1t::Stage1Layer2TauAlgorithmImpPP::JetIsolation(int et, int ieta, int iphi,
							const std::vector<l1t::Jet> & jets) const {

  for(JetBxCollection::const_iterator jet = jets.begin();
      jet != jets.end(); jet++) {

    if (ieta==jet->hwEta() && iphi==jet->hwPhi()){

      double isolation = (double) (jet->hwPt() - et);
      return isolation/et;
    }
  }

  // set output
  return 999.;
}


//  Find if the neighbor with the highest Et is N, E, S, or W
string l1t::Stage1Layer2TauAlgorithmImpPP::findNESW(int ieta, int iphi, int neta, int nphi) const {

  int deltaPhi = iphi - nphi;
  if (std::abs(deltaPhi) == L1CaloRegionDetId::N_PHI-1)
    deltaPhi = -deltaPhi/std::abs(deltaPhi); //18 regions in phi
//  std::cout << ieta << " " << neta << std::endl;
  int deltaEta = ieta - neta;
//  std::cout << "deltaPhi, deltaEta " << deltaPhi << ", " << deltaEta << std::endl;
  if ((std::abs(deltaPhi) +  std::abs(deltaEta)) < 2) {
    if (deltaEta==-1) {
      return "isEast";
    }
    else if (deltaEta==0) {
      if (deltaPhi==-1) {
	return "isNorth";
      }
      if (deltaPhi==1) {
	return "isSouth";
      }
      if (deltaPhi==0){
        return "isOverlap";
      }
    }
    else {
      return "isWest";
    }
  }

  return "999";

}
