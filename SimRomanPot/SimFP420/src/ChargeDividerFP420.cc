///////////////////////////////////////////////////////////////////////////////
// File: ChargeDividerFP420
// Date: 12.2006
// Description: ChargeDividerFP420  for FP420
// Modifications: std::
///////////////////////////////////////////////////////////////////////////////
#include "SimRomanPot/SimFP420/interface/ChargeDividerFP420.h"

using namespace std;
#include<vector>
//#define mydigidebug3



//DigitizerFP420::DigitizerFP420(const edm::ParameterSet& conf):conf_(conf),stripDigitizer_(new FP420DigiMain(conf)) 
ChargeDividerFP420::ChargeDividerFP420(double pit, double az420, double azD2, double azD3){
  //                           pit - is really moduleThickness here !!!
#ifdef mydigidebug3
cout << "ChargeDividerFP420.h: constructor" << endl;
cout << "peakMode = " << peakMode << "fluctuateCharge=   "<< fluctuateCharge <<  "chargedivisionsPerHit = "  << chargedivisionsPerHit << "deltaCut=   "<< deltaCut << endl;
#endif
  // Initialization:
  theFP420NumberingScheme = new FP420NumberingScheme();


  // Run APV in peak instead of deconvolution mode, which degrades the time resolution
//  peakMode=true ; //     APVpeakmode
  peakMode=false; //     APVconvolutionmode
  
  // Enable interstrip Landau fluctuations within a cluster.
  fluctuateCharge=true;   
  
  // Number of segments per strip into which charge is divided during simulation.
  // If large the precision of simulation improves.
  chargedivisionsPerHit=10; // = or =20
 
  // delta cutoff in MeV, has to be same as in OSCAR (0.120425 MeV corresponding // to 100um range for electrons)
  //SimpleConfigurable<double>  ChargeDividerFP420::deltaCut(0.120425,
  deltaCut=0.120425;  //  DeltaProductionCut

  pitchcur= pit;// pitchcur - is really moduleThickness here !!!

  // but position before Stations:
  z420 = az420;  // dist between centers of 1st and 2nd stations
  zD2 = azD2;  // dist between centers of 1st and 2nd stations
  zD3 = azD3;  // dist between centers of 1st and 3rd stations
        


  //                                                                                                                           .
  //                                                                                                                           .
  //  -300     -209.2             -150              -90.8                        0                                           +300
  //                                                                                                                           .
  //            X  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | X                        station                                          .
  //                   8*13.3+ 2*6 = 118.4                                    center                                           .
  //                                                                                                                           .
  zStationBegPos[0] = -150. - (118.4+10.)/2 + z420; // 10. -arbitrary
  zStationBegPos[1] = zStationBegPos[0]+zD2;
  zStationBegPos[2] = zStationBegPos[0]+zD3;
  zStationBegPos[3] = zStationBegPos[0]+2*zD3;

}

// Virtual destructor needed.
ChargeDividerFP420::~ChargeDividerFP420() { 
//  if(verbosity>0) {
//    std::cout << "Destroying a ChargeDividerFP420" << std::endl;
//  }
  delete theFP420NumberingScheme;
  
}  
CDividerFP420::ionization_type 
ChargeDividerFP420::divide(
			      const FP420G4Hit& hit, const double& pitchcur) {
  // !!!
  //                                       pitchcur - is really moduleThickness here !!!
  // !!!



  // sign "-" mean not the same as "+" for middle point !!!
  G4ThreeVector direction = hit.getExitLocalP() - hit.getEntryLocalP();  
  // direction.mag() - length or (size of path) of the hit; direction/direction.mag() - cosines of direction

#ifdef mydigidebug3
std::cout << " CDividerFP420::ChargeDividerFP420:divide: direction= " << direction << std::endl;
std::cout << " CDividerFP420::ChargeDividerFP420:divide: direction.mag = " << direction.mag() << std::endl;
std::cout << " obtained as ExitLocalP = " << hit.getExitLocalP() << "  -  "<<  "  EntryLocalP = "  << hit.getEntryLocalP() << std::endl;
std::cout << "  pitchcur= " << pitchcur << std::endl;
std::cout << "  peakMode = " << peakMode << "  fluctuateCharge=   "<< fluctuateCharge <<  "  chargedivisionsPerHit = "  << chargedivisionsPerHit << "  deltaCut=   "<< deltaCut << std::endl;
#endif

  int NumberOfSegmentation =  

//     (int)(1+chargedivisionsPerHit*fabs(direction.x())/pitchcur); // equidistant in X
//    (int)(1+chargedivisionsPerHit*fabs(direction.z())/pitchcur); // equidistant in Z, but why?

    (int)(1+chargedivisionsPerHit*direction.mag()/pitchcur); // equidistant over hit path

#ifdef mydigidebug3
std::cout << "NumberOfSegmentation= " << NumberOfSegmentation << std::endl;
#endif
 
//  float eLoss = hit.energyLoss();  // Eloss in GeV
  float eLoss = hit.getEnergyLoss();  // Eloss in GeV
#ifdef mydigidebug3
std::cout << "CDividerFP420::ChargeDividerFP420:divide: eLoss=  " << eLoss << std::endl;
#endif

  //
  //     return the energyLoss weighted CR-RC shape peaked at t0.(PeakShape)
  //     return the energyLoss weighted with a gaussian centered at t0 (DeconvolutionShape)
  float decSignal = TimeResponse(hit);
#ifdef mydigidebug3
std::cout << "CDividerFP420::ChargeDividerFP420:divide: decSignal=  " << decSignal << std::endl;
#endif
 
  ionization_type _ionization_points;

  _ionization_points.resize(NumberOfSegmentation);

  float energy;

  // Fluctuate charge in track subsegments
  float* eLossVector = new float[NumberOfSegmentation];
 
#ifdef mydigidebug3
std::cout << "CDividerFP420::ChargeDividerFP420:divide:  resize done; then, fluctuateCharge ? = " << fluctuateCharge << std::endl;
#endif
  if( fluctuateCharge ) {
    int pid = hit.getParticleType();
    float momentum = hit.getPabs();
    float length = direction.mag();  // length or (size of path) of the hit;
#ifdef mydigidebug3
std::cout << "pid= " << pid << "momentum= " << momentum << "eLoss= " << eLoss << "length= " << length << std::endl;
#endif
    fluctuateEloss(pid, momentum, eLoss, length, NumberOfSegmentation, eLossVector);   
  }
 
  for ( int i = 0; i != NumberOfSegmentation; i++) {
    if( fluctuateCharge ) {
      energy=eLossVector[i]*decSignal/eLoss;
      //      EnergySegmentFP420 edu(energy,hit.entryPoint()+float((i+0.5)/NumberOfSegmentation)*direction);//take energy value from vector eLossVector  
      EnergySegmentFP420 edu(energy,hit.getEntryLocalP()+float((i+0.5)/NumberOfSegmentation)*direction);//take energy value from vector eLossVector  
      _ionization_points[i] = edu; //save
    }else{
      energy=decSignal/float(NumberOfSegmentation);
      //      EnergySegmentFP420 edu(energy,hit.entryPoint()+float((i+0.5)/NumberOfSegmentation)*direction);//take energy value from eLoss average over n.segments 
      EnergySegmentFP420 edu(energy,hit.getEntryLocalP()+float((i+0.5)/NumberOfSegmentation)*direction);//take energy value from eLoss average over n.segments 
      _ionization_points[i] = edu; //save
    }
  }

#ifdef mydigidebug3
std::cout << "CDividerFP420::ChargeDividerFP420:divide:  !!!  RESULT !!!" << std::endl;
    std::cout <<  " _ionization_points size = " << _ionization_points.size() << std::endl;
  for(unsigned int i = 0; i < _ionization_points.size(); i++ ) {
    std::cout <<  " eLossVector[i] i = " << i << eLossVector[i] << std::endl;
  }
#endif
 
  delete[] eLossVector;
  return _ionization_points;
}
    
void ChargeDividerFP420::fluctuateEloss(int pid, float particleMomentum, 
				      float eloss, float length, 
				      int NumberOfSegs,float elossVector[]) {

  // Get dedx for this track
  float dedx;
#ifdef mydigidebug3
std::cout << "fluctuateEloss: eloss=  " << eloss << "length=  " << length << "NumberOfSegs=  " << NumberOfSegs << std::endl;
#endif
  if( length > 0.) dedx = eloss/length;
  else dedx = eloss;

  //  double particleMass = 139.57; // Mass in MeV, Assume pion
  double particleMass = 938.; // Mass in MeV, Assume proton   ----  AZ
  //  if( particleTable->getParticleData(pid) ) {  // Get mass from the PDTable
  //    particleMass = 1000. * particleTable->getParticleData(pid)->mass(); //Conv. GeV to MeV
  //  }

  float segmentLength = length/NumberOfSegs;

  // Generate charge fluctuations.
  float de=0.;
  float sum=0.;
  double segmentEloss = (1000.*eloss)/NumberOfSegs; //eloss in MeV
#ifdef mydigidebug3
std::cout << "segmentLength=  " << segmentLength << "segmentEloss=  " << segmentEloss << std::endl;
#endif

  for (int i=0;i<NumberOfSegs;i++) {
    // The G4 routine needs momentum in MeV, mass in Mev, delta-cut in MeV,
    // track segment length in mm(!!!), segment eloss in MeV 
    // Returns fluctuated eloss in MeV
    //    double deltaCutoff = deltaCut.value(); // the cutoff is sometimes redefined inside, so fix it.
    double deltaCutoff = deltaCut;
    de = fluctuate.SampleFluctuations(double(particleMomentum*1000.),
				      particleMass, deltaCutoff, 
				      double(segmentLength),
				      segmentEloss )/1000.; //convert to GeV
    elossVector[i]=de;
    sum +=de;
  }

#ifdef mydigidebug3
std::cout << "sum=  " << sum << std::endl;
#endif
  if(sum>0.) {  // If fluctuations give eloss>0.
    // Rescale to the same total eloss
    float ratio = eloss/sum;
    for (int ii=0;ii<NumberOfSegs;ii++) elossVector[ii]= ratio*elossVector[ii];
  } else {  // If fluctuations gives 0 eloss
    float averageEloss = eloss/NumberOfSegs;
    for (int ii=0;ii<NumberOfSegs;ii++) elossVector[ii]= averageEloss; 
  }
  return;
}

float ChargeDividerFP420::TimeResponse( const FP420G4Hit& hit ) {
  if (peakMode) {
#ifdef mydigidebug3
std::cout << "ChargeDividerFP420:TimeResponse: call of PeakShape" << std::endl;
#endif
    return this->PeakShape( hit );
  } else {
#ifdef mydigidebug3
std::cout << "ChargeDividerFP420:TimeResponse: call of DeconvolutionShape" << std::endl;
#endif
    return this->DeconvolutionShape( hit );
  }
}
float ChargeDividerFP420::PeakShape(const FP420G4Hit& hit){
  // 
  // Aim:     return the energyLoss weighted CR-RC shape peaked at t0.
  // 

    float xEntry = hit.getX() - hit.getVx();
    float yEntry = hit.getY() - hit.getVy();
    float zEntry = hit.getZ() - hit.getVz();

    unsigned int unitID = hit.getUnitID();
    //      int sScale = 20;
      int det, zside, sector, zmodule;
      FP420NumberingScheme::unpackFP420Index(unitID, det, zside, sector, zmodule);
// intindex is a continues numbering of FP420
//	  int zScale=2;  unsigned int intindex = sScale*(sector - 1)+zScale*(zmodule - 1)+zside;
//	 int zScale=10;	unsigned int intindex = sScale*(sector - 1)+zScale*(zside - 1)+zmodule;

  float RRR      = sqrt(xEntry*xEntry + yEntry*yEntry + zEntry*zEntry);
  float costheta = zEntry / RRR ;
  //  float theta    = acos(min(max(costheta,float(-1.)),float(1.)));
  //  float dist = hit.det().position().mag();
  //  float dist = hit.localPosition().mag();//AZ
  //  float dist = hit.getEntry().mag();
  //  float dist = hit.getEntryLocalP().mag();
  float dist     = (zStationBegPos[sector-1] - hit.getVz()) / costheta;
  dist     = dist/10.;// mm  --> cm as light velocity = 30 cm/ns
#ifdef mydigidebug3
std::cout << "sector=" << sector << std::endl;
std::cout << "zmodule=" << zmodule << std::endl;
std::cout << "zStationBegPos[sector-1]=" << zStationBegPos[sector-1] << std::endl;
std::cout << "RRR=" << RRR << std::endl;
std::cout << "costheta=" << costheta << std::endl;
std::cout << "unitID=" << unitID << std::endl;
//std::cout << "thetaEntry=" << thetaEntry << std::endl;
//std::cout << "my theta=" << theta*180./3.1415927 << std::endl;
std::cout << "dist found =" << dist << std::endl;
#endif

  // Time when read out relative to time hit produced.
  float t0 = dist/30.;  // light velocity = 30 cm/ns
  float SigmaShape = 52.17; 
  float tofNorm = (hit.getTof() - t0)/SigmaShape;

  float readTimeNorm = -tofNorm;
  // return the energyLoss weighted with CR-RC shape peaked at t0.
#ifdef mydigidebug3
std::cout << "ChargeDividerFP420:PeakShape::dist=" << dist << std::endl;
std::cout << "t0=" <<t0  << std::endl;
std::cout << "hit.getTof()=" << hit.getTof() << std::endl;
std::cout << "tofNorm=" << tofNorm << std::endl;
std::cout << "1 + readTimeNorm=" << 1 + readTimeNorm << std::endl;
std::cout << "hit.getEnergyLoss()=" << hit.getEnergyLoss()  << std::endl;
std::cout << "(1 + readTimeNorm)*exp(-readTimeNorm)=" << (1 + readTimeNorm)*exp(-readTimeNorm)  << std::endl;
std::cout << "return=" << hit.getEnergyLoss()*(1 + readTimeNorm)*exp(-readTimeNorm)  << std::endl;
#endif
  if (1 + readTimeNorm > 0) {
    //    return hit.energyLoss()*(1 + readTimeNorm)*exp(-readTimeNorm);
    return hit.getEnergyLoss()*(1 + readTimeNorm)*exp(-readTimeNorm);
  } else {
    return 0.;
  }
}

float ChargeDividerFP420::DeconvolutionShape(const FP420G4Hit& hit){
  // 
  // Aim:     return the energyLoss weighted with a gaussian centered at t0 
  // 
//    float thetaEntry = hit.getThetaAtEntry();   
//    float phiEntry = hit->getPhiAtEntry();

    float xEntry = hit.getX() - hit.getVx();
    float yEntry = hit.getY() - hit.getVy();
    float zEntry = hit.getZ() - hit.getVz();

    unsigned int unitID = hit.getUnitID();
    //      int sScale = 20;
      int det, zside, sector, zmodule;
      FP420NumberingScheme::unpackFP420Index(unitID, det, zside, sector, zmodule);
// intindex is a continues numbering of FP420
//	  int zScale=2;  unsigned int intindex = sScale*(sector - 1)+zScale*(zmodule - 1)+zside;
//	 int zScale=10;	unsigned int intindex = sScale*(sector - 1)+zScale*(zside - 1)+zmodule;

  float RRR      = sqrt(xEntry*xEntry + yEntry*yEntry + zEntry*zEntry);
  float costheta = zEntry / RRR ;
  //  float theta    = acos(min(max(costheta,float(-1.)),float(1.)));
  float dist     = (zStationBegPos[sector-1] - hit.getVz()) / costheta;
  dist     = dist/10.;// mm  --> cm as light velocity = 30 cm/ns
#ifdef mydigidebug3
std::cout << "sector=" << sector << std::endl;
std::cout << "zmodule=" << zmodule << std::endl;
std::cout << "zStationBegPos[sector-1]=" << zStationBegPos[sector-1] << std::endl;
std::cout << "RRR=" << RRR << std::endl;
std::cout << "costheta=" << costheta << std::endl;
std::cout << "unitID=" << unitID << std::endl;
//std::cout << "thetaEntry=" << thetaEntry << std::endl;
//std::cout << "my theta=" << theta*180./3.1415927 << std::endl;
std::cout << "dist found =" << dist << std::endl;
#endif

    float t0 = dist/30.;  // light velocity = 30 cm/ns
  float SigmaShape = 12.; 
//fun/pl 1*exp(-0.5*((0.1/30-x)/0.1)**2) 0. 0.08
//  float SigmaShape = 22.; 
  //  float tofNorm = (hit.tof() - t0)/SigmaShape;
  float tofNorm = (hit.getTof() - t0)/SigmaShape;
  // Time when read out relative to time hit produced.
  float readTimeNorm = -tofNorm;
  // return the energyLoss weighted with a gaussian centered at t0 
//  return hit.energyLoss()*exp(-0.5*readTimeNorm*readTimeNorm);
#ifdef mydigidebug3
std::cout << "ChargeDividerFP420:DeconvolutionShape::dist=" << dist << std::endl;
std::cout << "t0=" <<t0  << std::endl;
std::cout << "hit.getTof()=" << hit.getTof() << std::endl;
std::cout << "tofNorm=" << tofNorm << std::endl;
std::cout << "hit.getEnergyLoss()=" << hit.getEnergyLoss()  << std::endl;
std::cout << "exp(-0.5*readTimeNorm*readTimeNorm)=" << exp(-0.5*readTimeNorm*readTimeNorm)  << std::endl;
std::cout << "return=" << hit.getEnergyLoss()*exp(-0.5*readTimeNorm*readTimeNorm)  << std::endl;
#endif
  return hit.getEnergyLoss()*exp(-0.5*readTimeNorm*readTimeNorm);
}

