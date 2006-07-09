#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include <stdexcept>

ESDetId EcalPreshowerTopology::incrementIy(const ESDetId& id) const {
  try 
    {
      if (!(*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(id))
	{
	  return ESDetId(0);
	}
      ESDetId nextPoint;
      //Strips orientend along x direction for plane 2
      if (id.plane() == 2)
	{
	  if (id.strip() < 32 )
	    //Incrementing just strip number
	    nextPoint=ESDetId(id.strip()+1,id.six(),id.siy(),id.plane(),id.zside());
	  else
	    //Changing wafer
	    nextPoint=ESDetId(1,id.six(),id.siy()+1,id.plane(),id.zside());
	}
      //Strips orientend along y direction for plane 1
      else if (id.plane() == 1)
	{
	  //Changing wafer
	  nextPoint=ESDetId(id.strip(),id.six(),id.siy()+1,id.plane(),id.zside());
	}
      else
	return ESDetId(0);
      
      if ((*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(nextPoint))
	return nextPoint;
      else
	return ESDetId(0);
    } 
  catch ( std::runtime_error &e ) 
    { 
      //      std::cout << "    Catching: runtime_error " << e.what() << std::endl;
      return ESDetId(0);
    }
}


ESDetId EcalPreshowerTopology::decrementIy(const ESDetId& id) const {
  try 
    {
      if (!(*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(id))
	{
	  return ESDetId(0);
	}
      ESDetId nextPoint;
      //Strips orientend along x direction for plane 2
      if (id.plane() == 2)
	{
	  if (id.strip() >1 )
	    //Decrementing just strip number
	    nextPoint=ESDetId(id.strip()-1,id.six(),id.siy(),id.plane(),id.zside());
	  else
	    //Changing wafer
	    nextPoint=ESDetId(32,id.six(),id.siy()-1,id.plane(),id.zside());
	}
      //Strips orientend along y direction for plane 1
      else if (id.plane() == 1)
	{
	  //Changing wafer
	  nextPoint=ESDetId(id.strip(),id.six(),id.siy()-1,id.plane(),id.zside());
	}
      else
	return ESDetId(0);

      if ((*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(nextPoint))
	return nextPoint;
      else
	return ESDetId(0);
    } 
  catch ( std::runtime_error &e ) 
    { 
      return ESDetId(0);
    }
}

ESDetId EcalPreshowerTopology::incrementIx(const ESDetId& id) const {
  try 
    {
      if (!(*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(id))
	{
	  return ESDetId(0);
	}      
      ESDetId nextPoint;

      //Strips orientend along x direction for plane 2
      if (id.plane() == 2)
	{
	  //Changing wafer
	  nextPoint=ESDetId(id.strip(),id.six()+1,id.siy(),id.plane(),id.zside());
	}
      //Strips orientend along y direction for plane 1
      else if (id.plane() == 1)
	{
	  if (id.strip() < 32 )
	    //Incrementing just strip number
	    nextPoint=ESDetId(id.strip()+1,id.six(),id.siy(),id.plane(),id.zside());
	  else
	    //Changing wafer
	    nextPoint=ESDetId(1,id.six()+1,id.siy(),id.plane(),id.zside());
	}
      else
	return ESDetId(0);

      if ((*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(nextPoint))
	return nextPoint;
      else
	return ESDetId(0);
    } 
  catch ( std::runtime_error &e ) 
    { 
      return ESDetId(0);
    }
}

ESDetId EcalPreshowerTopology::decrementIx(const ESDetId& id) const {
  try 
    {
      if (!(*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(id))
	{
	  return ESDetId(0);
	}      
      ESDetId nextPoint;
      //Strips orientend along x direction for plane 2
      if (id.plane() == 2)
	{
	  //Changing wafer
	  nextPoint=ESDetId(id.strip(),id.six()-1,id.siy(),id.plane(),id.zside());
	}
      //Strips orientend along y direction for plane 1
      else if (id.plane() == 1)
	{
	  if (id.strip() > 1 )
	    //Decrementing just strip number
	    nextPoint=ESDetId(id.strip()-1,id.six(),id.siy(),id.plane(),id.zside());
	  else
	    //Changing wafer
	    nextPoint=ESDetId(32,id.six()-1,id.siy(),id.plane(),id.zside());
	}
      else
	return ESDetId(0);

      if ((*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(nextPoint))
	return nextPoint;
      else
	return ESDetId(0);
    } 
  catch ( std::runtime_error &e ) 
    { 
      return ESDetId(0);
    }
}


ESDetId EcalPreshowerTopology::incrementIz(const ESDetId& id) const {
  try 
    {
      if (!(*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(id))
	{
	  return ESDetId(0);
	}      
      ESDetId nextPoint;
      nextPoint=ESDetId(id.strip(),id.six(),id.siy(),id.plane()+1,id.zside());
      if ((*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(nextPoint))
	return nextPoint;
      else
	return ESDetId(0);
    } 
  catch ( std::runtime_error &e ) 
    { 
      return ESDetId(0);
    }
}


ESDetId EcalPreshowerTopology::decrementIz(const ESDetId& id) const {
  try 
    {
      if (!(*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(id))
	{
	  return ESDetId(0);
	}      
      ESDetId nextPoint;
      nextPoint=ESDetId(id.strip(),id.six(),id.siy(),id.plane()-1,id.zside());
      if ((*theGeom_).getSubdetectorGeometry(DetId::Ecal,EcalPreshower)->present(nextPoint))
	return nextPoint;
      else
	return ESDetId(0);
    } 
  catch ( std::runtime_error &e ) 
    { 
      return ESDetId(0);
    }
}






