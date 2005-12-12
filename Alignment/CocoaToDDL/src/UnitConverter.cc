using namespace std;

#include "OpticalAlignment/CocoaToDDL/interface/UnitConverter.h"
#include "OpticalAlignment/CocoaToDDL/interface/CocoaUnitsTable.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include <sstream>
#include <strstream>
#include <string>

/*
UnitConverter::UnitConverter(const G4BestUnit & bu)
 : bu_(bu)
{
  //ostrstream s;
  //s << bu;
  
}
*/

UnitConverter::UnitConverter(ALIdouble val, const ALIstring & category)
 : bu_(new CocoaBestUnit(val,category)), angl_(false)
{ 
   if (category=="Angle")
    angl_=true;
}

UnitConverter::~UnitConverter()
{
  delete bu_;
}  

std::string UnitConverter::ucstring()
{

   std::ostrstream str;
   
   if (angl_) {
     str.precision(11);
     double x = ( *(bu_->GetValue()) ) / deg;
     str << x << string("*deg") << '\0';
     return string(str.str());
   
   }
   else {
     str << *bu_ << '\0';
     std::string s(str.str());
     return s.replace(s.find(" "),1,"*");
   }  
   //return s;

}

/*
ostream & operator<<(ostream & os, const UnitConverter & uc)
{
  ostrstream temp;
  //temp << uc.bu_;
  //temp << '\0';
  //string s(temp.str());
  //cout << "NOW: " << s << endl;
  os << *(uc.bu_);
}

*/
