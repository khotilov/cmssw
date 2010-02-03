#include "DetectorDescription/Core/interface/DDCompactViewImpl.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

// Message logger.
#include "FWCore/MessageLogger/interface/MessageLogger.h"

DDCompactViewImpl::DDCompactViewImpl(const DDLogicalPart & rootnodedata)
  :  root_(rootnodedata)//, rootWalker_(0)
{
  //  edm::LogInfo("DDCompactViewImpl") << "Root node data = " << rootnodedata << std::endl;
   //buildTree(rootnodedata,root_);
   //buildGraph();
}


DDCompactViewImpl::DDCompactViewImpl()
 : root_("")
{

}


DDCompactViewImpl::~DDCompactViewImpl() 
{  
   GraphNav::adj_list::size_type it = 0;
   for (; it < graph_.size() ; ++it) {
     GraphNav::edge_range erange = graph_.edges(it); //it->second.begin();
     for(; erange.first != erange.second; ++(erange.first)) {
       DDPosData * pd = graph_.edgeData(erange.first->second);
       delete pd;
       pd=0;
     }  
   }
   edm::LogInfo("DDCompactViewImpl") << std::endl << "DDD transient representation has been destructed." << std::endl << std::endl;   
}


// deprecated
void DDCompactViewImpl::buildPaths()
{
  //paths_ = new GraphNavPaths(graph_,root_); // good luck!
}


/*
std::pair<bool,DDPhysicalPart> DDCompactViewImpl::goTo(const DDPartSelector & path) const
{
  std::pair<bool,DDPhysicalPart> result(std::make_pair(false,DDPhysicalPart()));
  DDPhysicalPart & pp(result.second);
  DDTranslation trans;
  DDRotationMatrix rot;
  return result;
}

*/

graphwalker<DDLogicalPart,DDPosData*> DDCompactViewImpl::walker() const
{
   DCOUT('C',"DDCompactView::walker() root_=" << root_);
   return graphwalker<DDLogicalPart,DDPosData*>(graph_,root_);
}

// calculates the weight and caches it in LogicalPartImpl
//double DDCompactViewImpl::weight(DDLogicalPart & part)
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDMaterial.h"
double DDCompactViewImpl::weight(const DDLogicalPart & aPart) const
{
 // return 0;

   if (!aPart)
     return -1;
   DDLogicalPart part = aPart;
   //double result;  
   if (part.weight())
     return part.weight();
   
   // weigth = (parent.vol - children.vol)*parent.density + weight(children)
   double childrenVols=0;
   double childrenWeights=0;
   walker_type walker(graph_,part);
   if(walker.firstChild()) {
     bool doIt=true;
     while(doIt) {
       double a_vol;
       DDLogicalPart child(walker.current().first);
       a_vol=child.solid().volume();
       if (a_vol <= 0.) {
         edm::LogError("DDCompactViewImpl")  << "DD-WARNING: volume of solid=" << aPart.solid() 
	       << "is zero or negative, vol=" << a_vol/m3 << "m3" << std::endl;
       }
       DCOUT_V('C', "DC: weightcalc, currently=" << child.ddname().name()
            << " vol=" << a_vol/cm3 << "cm3");
       childrenVols += a_vol;
       childrenWeights += weight(child); // recursive
       doIt=walker.nextSibling();
     }
   }
   
   double dens = part.material().density();
   if (dens <=0) {
     edm::LogError("DDCompactViewImpl")  << "DD-WARNING: density of material=" << part.material().ddname() 
	   << " is negative or zero, rho=" << dens/g*cm3 << "g/cm3" << std::endl;
   }
   double p_vol  = part.solid().volume();
   double w =   (p_vol - childrenVols)*dens + childrenWeights; 
   if (  (fabs(p_vol) - fabs(childrenVols))/fabs(p_vol) > 1.01 ) {
     edm::LogError("DDCompactViewImpl")  << "DD-WARNING: parent-volume smaller than children, parent=" 
          << part.ddname() << " difference-vol=" 
	   << (p_vol - childrenVols)/m3 << "m3, this is " 
	   << (childrenVols - p_vol)/p_vol << "% of the parent-vol." << std::endl;
   }
  
   //part.rep().weight_=w;
   part.weight() = w;
   return w;
   
}

void DDCompactViewImpl::swap( DDCompactViewImpl& implToSwap ) {
  //  root_.rep().swap(implToSwap.root_.rep()); // I do not want to do this because to construct a DDCompactView that is NOT global one needs to do DDCompactView(lp)
  graph_.swap(implToSwap.graph_);
}

/*
expnode_t * DDCompactViewImpl::expand(const DDPartSelector & path) const
{
   //FIXME: DDCompactViewImpl::expand - add?? isUnique(path) ??...
   expnode_t * result = 0;
   return result;
    
}

*/
