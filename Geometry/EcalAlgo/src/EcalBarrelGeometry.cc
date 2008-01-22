#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include <CLHEP/Geometry/Point3D.h>
#include <CLHEP/Geometry/Plane3D.h>

#include <iomanip>

using namespace std;

EcalBarrelGeometry::EcalBarrelGeometry() :_nnxtalEta(0),
					  _nnxtalPhi(0 ),
					  _PhiBaskets(0)
{
}


EcalBarrelGeometry::~EcalBarrelGeometry() 
{
}

// Get closest cell, etc...
DetId 
EcalBarrelGeometry::getClosestCell(const GlobalPoint& r) const 
{

  // z is the easy one
  int leverx = 1;
  int levery = 1;
  double pointz = r.z();
  int zbin=1;
  if(pointz<0)
    zbin=-1;

  // Now find the closest eta
  double pointeta = r.eta();
  //  double eta;
  double deta=999.;
  int etabin=1;
  
  int guessed_eta = (int)( fabs(pointeta) / 0.0174)+1;
  int guessed_eta_begin = guessed_eta-1;
  int guessed_eta_end   = guessed_eta+1;
  if (guessed_eta_begin < 1) guessed_eta_begin = 1;
  if (guessed_eta_end > 85) guessed_eta_end = 85;

    for(int bin=guessed_eta_begin; bin<= guessed_eta_end; bin++)
      {
	try
	  {
	    if (!present(EBDetId(zbin*bin,1,EBDetId::ETAPHIMODE)))
	      continue;

	    double eta = getGeometry(EBDetId(zbin*bin,1,EBDetId::ETAPHIMODE))->getPosition().eta();

	    if(fabs(pointeta-eta)<deta)
	      {
		deta=fabs(pointeta-eta);
		etabin=bin;
	      }
	    else break;
	  }
	catch ( cms::Exception &e ) 
	  {
	  }
      }
    

  // Now the closest phi. always same number of phi bins(!?)
  const double twopi = M_PI+M_PI;

  // 10 degree tilt
  const double tilt=twopi/36.;
  double pointphi = r.phi()+tilt;

  // put phi in correct range (0->2pi)
  if(pointphi > twopi)
    pointphi -= twopi;
  if(pointphi < 0)
    pointphi += twopi;

  //calculate phi bin, distinguish + and - eta
  int phibin = static_cast<int>(pointphi / (twopi/_nnxtalPhi)) + 1;
  //   if(point.z()<0.0)
  //     {
  //       phibin = nxtalPhi/2 - 1 - phibin;
  //       if(phibin<0)
  //         phibin += nxtalPhi;
  //     }
  try
    {
      EBDetId myCell(zbin*etabin,phibin,EBDetId::ETAPHIMODE);

      if (!present(myCell))
	return DetId(0);
      
      HepPoint3D  A;
      HepPoint3D  B;
      HepPoint3D  C;
      HepPoint3D  point(r.x(),r.y(),r.z());

      // D.K. : equation of plane : AA*x+BB*y+CC*z+DD=0;
      // finding equation for each edge

      // Since the point can lie between crystals, it is necessary to keep track of the movements
      // to avoid infinite loops
      std::vector<double> history;
      history.resize(4,0.);
      //
      // stop movement in eta direction when closest cell was found (point between crystals)
      int start = 1;
      int counter = 0;
      // Moving until find closest crystal in eta and phi directions (leverx and levery)
      while  (leverx==1 || levery == 1)
	{
	  leverx = 0;
	  levery = 0;
	  const CaloCellGeometry::CornersVec& corners 
	     ( getGeometry(myCell)->getCorners() ) ;
	  std::vector<double> SS;

	  // compute the distance of the point with respect of the 4 crystal lateral planes
	  for (short i=0; i < 4 ; ++i)
	    {
	      A = HepPoint3D(corners[i%4].x(),corners[i%4].y(),corners[i%4].z());
	      B = HepPoint3D(corners[(i+1)%4].x(),corners[(i+1)%4].y(),corners[(i+1)%4].z());
	      C = HepPoint3D(corners[4+(i+1)%4].x(),corners[4+(i+1)%4].y(),corners[4+(i+1)%4].z());
	      HepPlane3D plane(A,B,C);
	      plane.normalize();
	      double distance = plane.distance(point);
	      if(plane.d()>0.) distance=-distance;
	      if (corners[0].z()<0.) distance=-distance;
	      SS.push_back(distance);
	    }

	  // SS's - normals
	  // check position of the point with respect to opposite side of crystal
	  // if SS's have opposite sign, the  point lies inside that crystal

	  if ( ( SS[0]>0.&&SS[2]>0. )||( SS[0]<0.&&SS[2]<0. ) )
	    {
	      levery = 1;
	      if ( history[0]>0. && history[2]>0. && SS[0]<0 && SS[2]<0 &&
		   (abs(SS[0])+abs(SS[2]))> (abs(history[0])+abs(history[2]))) levery = 0  ;
	      if ( history[0]<0. && history[2]<0. && SS[0]>0 && SS[2]>0 &&
		   (abs(SS[0])+abs(SS[2]))> (abs(history[0])+abs(history[2]))) levery = 0  ;


	      if (SS[0]>0. )
		{
		  EBDetId nextPoint;
		  if (myCell.iphi()==EBDetId::MIN_IPHI) 
		    nextPoint=EBDetId(myCell.ieta(),EBDetId::MAX_IPHI);
		  else 
		    nextPoint=EBDetId(myCell.ieta(),myCell.iphi()-1);
		  if (present(nextPoint))
		    myCell=nextPoint;
		  else
		    levery=0;		  
		}
	      else
		{
		  EBDetId nextPoint;
		  if (myCell.iphi()==EBDetId::MAX_IPHI)
		    nextPoint=EBDetId(myCell.ieta(),EBDetId::MIN_IPHI);
		  else
		    nextPoint=EBDetId(myCell.ieta(),myCell.iphi()+1);
		  if (present(nextPoint))
		    myCell=nextPoint;
		  else
		    levery=0;
		}
	    }


	  if ( ( ( SS[1]>0.&&SS[3]>0. )||( SS[1]<0.&&SS[3]<0. )) && start==1  )
	    {
	      leverx = 1;

	      if ( history[1]>0. && history[3]>0. && SS[1]<0 && SS[3]<0 &&
		   (abs(SS[1])+abs(SS[3]))> (abs(history[1])+abs(history[3])) )
		{
		  leverx = 0;
		  start = 0;
		}

	      if ( history[1]<0. && history[3]<0. && SS[1]>0 && SS[3]>0 &&
		   (abs(SS[1])+abs(SS[3]))> (abs(history[1])+abs(history[3])) )
		{
		  leverx = 0;
		  start = 0;
		}


	      if (SS[1]>0.)
		{
		  EBDetId nextPoint;
		  if (myCell.ieta()==-1) 
		    nextPoint=EBDetId (1,myCell.iphi());
		  else 
		    {
		      int nieta= myCell.ieta()+1;
		      if(nieta==86) nieta=85;
		      nextPoint=EBDetId(nieta,myCell.iphi());
		    }
		  if (present(nextPoint))
		    myCell = nextPoint;
		  else
		    leverx = 0;
		}
	      else
		{
		  EBDetId nextPoint;
		  if (myCell.ieta()==1) 
		    nextPoint=EBDetId(-1,myCell.iphi());
		  else 
		    {
		      int nieta=myCell.ieta()-1;
		      if(nieta==-86) nieta=-85;
		      nextPoint=EBDetId(nieta,myCell.iphi());
		    }
		  if (present(nextPoint))
		    myCell = nextPoint;
		  else
		    leverx = 0;
		}
	    }
	  
	  // Update the history. If the point lies between crystals, the closest one
	  // is returned
	  history =SS;
	  
	  counter++;
	  if (counter == 10)
	    {
	      leverx=0;
	      levery=0;
	    }
	}
      // D.K. if point lies netween cells, take a closest cell.
      return DetId(myCell);
    }
  catch ( cms::Exception &e ) 
    { 
      return DetId(0);
    }

}

#include <iostream>
CaloSubdetectorGeometry::DetIdSet EcalBarrelGeometry::getCells(const GlobalPoint& r, double dR) const {
  double lowEta=r.eta()-dR;
  double highEta=r.eta()+dR;
  
  CaloSubdetectorGeometry::DetIdSet dis;
  if (highEta<-2.0 || lowEta>2.0) return dis;

  const double scale=EBDetId::MAX_IPHI/(2*M_PI);
  int ieta_center=int(r.eta()*scale+((r.z()<0)?(-1):(1)));
  double phi=r.phi();  if (phi<0) phi+=2*M_PI;
  int iphi_center=int(phi*scale+0.5)+10; // -10 for current geometry!
  if (iphi_center<=0) iphi_center+=EBDetId::MAX_IPHI;
  if (iphi_center>EBDetId::MAX_IPHI) iphi_center-=EBDetId::MAX_IPHI;

  double fr=dR/(2*M_PI/360);
  int idr=int(fr+0.5)+1;
  int idr2=int(fr*fr+4*fr+0.5);
  int idr2m=int(fr*fr-4*fr+0.5);
  for (int de=-idr; de<=idr; de++) {
    int ieta=de+ieta_center;

    if (ieta<=0 && ieta_center>0) ieta-=1; // go a little further
    if (ieta>=0 && ieta_center<0) ieta+=1; // go a little further
	
    if (ieta<-EBDetId::MAX_IETA || ieta==0 || ieta>EBDetId::MAX_IETA) continue; // not in EB
    for (int dp=-idr; dp<=idr; dp++) {
      int irange2=dp*dp+de*de;
	  
      int iphi=iphi_center+dp;
      if (iphi<1) iphi+=EBDetId::MAX_IPHI;
      else if (iphi>EBDetId::MAX_IPHI) iphi-=EBDetId::MAX_IPHI;
      
      if (irange2>idr2) continue;

      EBDetId id(ieta,iphi);
      bool ok=(irange2<idr2m);
      if (!ok) {
	const CaloCellGeometry* cell=getGeometry(id);
	ok=(cell!=0 && deltaR(r,cell->getPosition())<=dR);
	if (cell==0) std::cout << id << std::endl;	  
      }
      if (ok) dis.insert(id);
    }
  }
  return dis;
}
