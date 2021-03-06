// -*- C++ -*-
//
// Package:    DeDxEstimatorProducer
// Class:      DeDxEstimatorProducer
// 
/**\class DeDxEstimatorProducer DeDxEstimatorProducer.cc RecoTracker/DeDxEstimatorProducer/src/DeDxEstimatorProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  andrea
//         Created:  Thu May 31 14:09:02 CEST 2007
//    Code Updates:  loic Quertenmont (querten)
//         Created:  Thu May 10 14:09:02 CEST 2008
//
//


// system include files
#include "RecoTracker/DeDx/plugins/DeDxEstimatorProducer.h"


using namespace reco;
using namespace std;
using namespace edm;


void DeDxEstimatorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<string>("estimator","generic");
  desc.add<edm::InputTag>("tracks",edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("trajectoryTrackAssociation",edm::InputTag("generalTracks"));
  desc.add<bool>("UseTrajectory",true);  
  desc.add<bool>("UsePixel",false); 
  desc.add<bool>("UseStrip",true); 
  desc.add<double>("MeVperADCPixel",3.61e-06*265);
  desc.add<double>("MeVperADCStrip",3.61e-06);
  desc.add<bool>("ShapeTest",true);      
  desc.add<bool>("UseCalibration",false);  
  desc.add<string>("calibrationPath", "");
  desc.add<string>("Reccord", "SiStripDeDxMip_3D_Rcd");
  desc.add<string>("ProbabilityMode", "Accumulation");
  desc.add<double>("fraction", 0.4);
  desc.add<double>("exponent",-2.0);

  descriptions.add("DeDxEstimatorProducer",desc);
}


DeDxEstimatorProducer::DeDxEstimatorProducer(const edm::ParameterSet& iConfig)
{

   produces<ValueMap<DeDxData> >();

   string estimatorName = iConfig.getParameter<string>("estimator");
   if     (estimatorName == "median")              m_estimator = new MedianDeDxEstimator(iConfig);
   else if(estimatorName == "generic")             m_estimator = new GenericAverageDeDxEstimator  (iConfig);
   else if(estimatorName == "truncated")           m_estimator = new TruncatedAverageDeDxEstimator(iConfig);
   else if(estimatorName == "unbinnedFit")         m_estimator = new UnbinnedFitDeDxEstimator(iConfig);
   else if(estimatorName == "productDiscrim")      m_estimator = new ProductDeDxDiscriminator(iConfig);
   else if(estimatorName == "btagDiscrim")         m_estimator = new BTagLikeDeDxDiscriminator(iConfig);
   else if(estimatorName == "smirnovDiscrim")      m_estimator = new SmirnovDeDxDiscriminator(iConfig);
   else if(estimatorName == "asmirnovDiscrim")     m_estimator = new ASmirnovDeDxDiscriminator(iConfig);

  //Commented for now, might be used in the future
//   MaxNrStrips         = iConfig.getUntrackedParameter<unsigned>("maxNrStrips"        ,  255);

   m_tracksTag = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
   m_trajTrackAssociationTag   = consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("trajectoryTrackAssociation"));
   useTrajectory = iConfig.getParameter<bool>("UseTrajectory"); 

   usePixel = iConfig.getParameter<bool>("UsePixel"); 
   useStrip = iConfig.getParameter<bool>("UseStrip");
   meVperADCPixel = iConfig.getParameter<double>("MeVperADCPixel"); 
   meVperADCStrip = iConfig.getParameter<double>("MeVperADCStrip"); 

   shapetest = iConfig.getParameter<bool>("ShapeTest");
   useCalibration = iConfig.getParameter<bool>("UseCalibration");
   m_calibrationPath = iConfig.getParameter<string>("calibrationPath");

   if(!usePixel && !useStrip)
   edm::LogWarning("DeDxHitsProducer") << "Pixel Hits AND Strip Hits will not be used to estimate dEdx --> BUG, Please Update the config file";

}


DeDxEstimatorProducer::~DeDxEstimatorProducer()
{
  delete m_estimator;
}

// ------------ method called once each job just before starting event loop  ------------
void  DeDxEstimatorProducer::beginRun(edm::Run const& run, const edm::EventSetup& iSetup)
{
   if(useCalibration && calibGains.size()==0){
      edm::ESHandle<TrackerGeometry> tkGeom;
      iSetup.get<TrackerDigiGeometryRecord>().get( tkGeom );
      m_off = tkGeom->offsetDU(GeomDetEnumerators::PixelBarrel); //index start at the first pixel

      DeDxTools::makeCalibrationMap(m_calibrationPath, *tkGeom, calibGains, m_off);
   }

   m_estimator->beginRun(run, iSetup);
}



void DeDxEstimatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  auto_ptr<ValueMap<DeDxData> > trackDeDxEstimateAssociation(new ValueMap<DeDxData> );  
  ValueMap<DeDxData>::Filler filler(*trackDeDxEstimateAssociation);

  edm::Handle<reco::TrackCollection> trackCollectionHandle;
  iEvent.getByToken(m_tracksTag,trackCollectionHandle);

  Handle<TrajTrackAssociationCollection> trajTrackAssociationHandle;
  if(useTrajectory)iEvent.getByToken(m_trajTrackAssociationTag, trajTrackAssociationHandle);

  std::vector<DeDxData> dedxEstimate( trackCollectionHandle->size() );

  TrajTrackAssociationCollection::const_iterator cit;
  if(useTrajectory)cit = trajTrackAssociationHandle->begin();
  for(unsigned int j=0;j<trackCollectionHandle->size();j++){            
     const reco::TrackRef track = reco::TrackRef( trackCollectionHandle.product(), j );

     int NClusterSaturating = 0; 
     DeDxHitCollection dedxHits;
   
     if(useTrajectory){  //trajectory allows to take into account the local direction of the particle on the module sensor --> muc much better 'dx' measurement
        const edm::Ref<std::vector<Trajectory> > traj = cit->key; cit++;
        const vector<TrajectoryMeasurement> & measurements = traj->measurements();
        for(vector<TrajectoryMeasurement>::const_iterator it = measurements.begin(); it!=measurements.end(); it++){
           TrajectoryStateOnSurface trajState=it->updatedState();
           if( !trajState.isValid()) continue;
     
           const TrackingRecHit * recHit=(*it->recHit()).hit();
           if(!recHit || !recHit->isValid())continue;
           LocalVector trackDirection = trajState.localDirection();
           float cosine = trackDirection.z()/trackDirection.mag();

           processHit(recHit, trajState.localMomentum().mag(), cosine, dedxHits, NClusterSaturating);
        }

     }else{ //assume that the particles trajectory is a straight line originating from the center of the detector  (can be improved)
        for(unsigned int h=0;h<track->recHitsSize();h++){
           const TrackingRecHit* recHit = &(*(track->recHit(h)));
           if(!recHit || !recHit->isValid())continue;
           auto const & thit = static_cast<BaseTrackerRecHit const&>(*recHit);
           if(!thit.isValid())continue;//make sure it's a tracker hit

           const GlobalVector& ModuleNormal = recHit->detUnit()->surface().normalVector();         
           float cosine = (track->px()*ModuleNormal.x()+track->py()*ModuleNormal.y()+track->pz()*ModuleNormal.z())/track->p();

           processHit(recHit, track->p(), cosine, dedxHits, NClusterSaturating);
        } 
     }

     sort(dedxHits.begin(),dedxHits.end(),less<DeDxHit>());   
     std::pair<float,float> val_and_error = m_estimator->dedx(dedxHits);

     //WARNING: Since the dEdX Error is not properly computed for the moment
     //It was decided to store the number of saturating cluster in that dataformat
     val_and_error.second = NClusterSaturating; 
     dedxEstimate[j] = DeDxData(val_and_error.first, val_and_error.second, dedxHits.size() );
  }
  ///////////////////////////////////////

  filler.insert(trackCollectionHandle, dedxEstimate.begin(), dedxEstimate.end());

  // fill the association map and put it into the event
  filler.fill();
  iEvent.put(trackDeDxEstimateAssociation);   
}


void DeDxEstimatorProducer::processHit(const TrackingRecHit* recHit, float trackMomentum, float& cosine, reco::DeDxHitCollection& dedxHits, int& NClusterSaturating){
      auto const & thit = static_cast<BaseTrackerRecHit const&>(*recHit);
      if(!thit.isValid())return;

      auto const & clus = thit.firstClusterRef();
      if(!clus.isValid())return;

      if(clus.isPixel()){
          if(!usePixel) return;

          auto& detUnit     = *(recHit->detUnit());
          float pathLen     = detUnit.surface().bounds().thickness()/fabs(cosine);
          float chargeAbs   = clus.pixelCluster().charge();
          float charge      = meVperADCPixel*chargeAbs/pathLen;
          dedxHits.push_back( DeDxHit( charge, trackMomentum, pathLen, thit.geographicalId()) );
       }else if(clus.isStrip() && !thit.isMatched()){
          if(!useStrip) return;

          auto& detUnit     = *(recHit->detUnit());
          int   NSaturating = 0;
          float pathLen     = detUnit.surface().bounds().thickness()/fabs(cosine);
          float chargeAbs   = DeDxTools::getCharge(&(clus.stripCluster()),NSaturating, detUnit, calibGains, m_off);
          float charge      = meVperADCStrip*chargeAbs/pathLen;
          if(!shapetest || (shapetest && DeDxTools::shapeSelection(clus.stripCluster().amplitudes()))){
             dedxHits.push_back( DeDxHit( charge, trackMomentum, pathLen, thit.geographicalId()) );
             if(NSaturating>0)NClusterSaturating++;
          }
       }else if(clus.isStrip() && thit.isMatched()){
          if(!useStrip) return;
          const SiStripMatchedRecHit2D* matchedHit=dynamic_cast<const SiStripMatchedRecHit2D*>(recHit);
          if(!matchedHit)return;
          const GluedGeomDet* gdet = static_cast<const GluedGeomDet*>(matchedHit->det());

          auto& detUnitM     = *(gdet->monoDet());
          int   NSaturating = 0;
          float pathLen     = detUnitM.surface().bounds().thickness()/fabs(cosine);
          float chargeAbs   = DeDxTools::getCharge(&(matchedHit->monoCluster()),NSaturating, detUnitM, calibGains, m_off);
          float charge      = meVperADCStrip*chargeAbs/pathLen;
          if(!shapetest || (shapetest && DeDxTools::shapeSelection(matchedHit->monoCluster().amplitudes()))){
             dedxHits.push_back( DeDxHit( charge, trackMomentum, pathLen, matchedHit->monoId()) );
             if(NSaturating>0)NClusterSaturating++;
          }

          auto& detUnitS     = *(gdet->stereoDet());
          NSaturating = 0;
          pathLen     = detUnitS.surface().bounds().thickness()/fabs(cosine);
          chargeAbs   = DeDxTools::getCharge(&(matchedHit->stereoCluster()),NSaturating, detUnitS, calibGains, m_off);
          charge      = meVperADCStrip*chargeAbs/pathLen;
          if(!shapetest || (shapetest && DeDxTools::shapeSelection(matchedHit->stereoCluster().amplitudes()))){
             dedxHits.push_back( DeDxHit( charge, trackMomentum, pathLen, matchedHit->stereoId()) );
             if(NSaturating>0)NClusterSaturating++;
          }      
       }
}



//define this as a plug-in
DEFINE_FWK_MODULE(DeDxEstimatorProducer);
