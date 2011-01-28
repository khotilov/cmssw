#include "FWPFEcalRecHitLegoProxyBuilder.h"

//______________________________________________________________________________________________________
void
FWPFEcalRecHitLegoProxyBuilder::scaleProduct( TEveElementList *parent, FWViewType::EType type, const FWViewContext *vc )
{
   FWViewEnergyScale *caloScale = vc->getEnergyScale();
   bool b = caloScale->getPlotEt();
   float maxVal = getMaxVal( b );
   typedef std::vector<FWPFLegoRecHit*> rh;

   // printf("FWPFEcalRecHitLegoProxyBuilder::scaleProduct >> scale %f \n", caloScale->getValToHeight());
   for( rh::iterator i = m_recHits.begin(); i != m_recHits.end(); ++i )
   {  // Tallest tower needs deciding still
      if( (*i)->isTallest() == false && (*i)->getEtEnergy( b ) == maxVal )
         (*i)->setIsTallest( true );

      (*i)->updateScale( vc);
   }
}

//______________________________________________________________________________________________________
void
FWPFEcalRecHitLegoProxyBuilder::localModelChanges( const FWModelId &iId, TEveElement *parent,
                                                   FWViewType::EType viewType, const FWViewContext *vc )
{
   for (TEveElement::List_i i = parent->BeginChildren(); i!= parent->EndChildren(); ++i)
   {
      {
         TEveStraightLineSet* line = dynamic_cast<TEveStraightLineSet*>(*i);
         if (line)
         {
            const FWDisplayProperties &p = item()->modelInfo( iId.index() ).displayProperties();
            line->SetMarkerColor( p.color() );
         }
      }
   }
}

//______________________________________________________________________________________________________
TEveVector
FWPFEcalRecHitLegoProxyBuilder::calculateCentre( const std::vector<TEveVector> &corners ) const
{
   TEveVector centre;

   for( size_t i = 0; i < corners.size(); ++i )
   {
      centre.fX += corners[i].fX;
      centre.fY += corners[i].fY;            // Get total for x,y,z values
      centre.fZ += corners[i].fZ;
   }
   centre *= 1.f / 8.f;

   return centre;   
}

//______________________________________________________________________________________________________
float
FWPFEcalRecHitLegoProxyBuilder::calculateEt( const TEveVector &centre, float E )
{
   TEveVector vec = centre;
   float et;

   vec.Normalize();
   vec *= E;
   et = vec.Perp();

   return et;
}

//______________________________________________________________________________________________________
void
FWPFEcalRecHitLegoProxyBuilder::build( const FWEventItem *iItem, TEveElementList *product, const FWViewContext *vc )
{
   float maxEnergy = 0.0f;
   float maxEt = 0.0f;
   size_t itemSize = iItem->size(); //cache size

   for( size_t index = 0; index < itemSize; ++index )
   {
      TEveCompound *itemHolder = createCompound();
      product->AddElement( itemHolder );

      const EcalRecHit &iData = modelData( index );
      const float *corners = item()->getGeom()->getCorners( iData.detid() );
      float energy, et;
      std::vector<TEveVector> etaphiCorners(8);
      TEveVector centre;

      if( corners == 0 )
         continue;

      int k = 3;
      for( int i = 0; i < 4; ++i )
      {
         int j = k * 3;
         TEveVector cv = TEveVector( corners[j], corners[j+1], corners[j+2] );
         etaphiCorners[i].fX = cv.Eta();                                     // Conversion of rechit X/Y values for plotting in Eta/Phi
         etaphiCorners[i].fY = cv.Phi();
         etaphiCorners[i].fZ = 0.0;

         etaphiCorners[i+4].fX = etaphiCorners[i].fX;                        // Top can simply be plotted exactly over the top of the bottom face
         etaphiCorners[i+4].fY = etaphiCorners[i].fY;
         etaphiCorners[i+4].fZ = 0.001;
         // printf("%f %f %d \n",  etaphiCorners[i].fX, etaphiCorners[i].fY, i);
         --k;
      }

      centre = calculateCentre( etaphiCorners );
      energy = iData.energy();
      et = calculateEt( centre, energy );
      context().voteMaxEtAndEnergy( et, energy );

      if( energy > maxEnergy )
         maxEnergy = energy;
      if( energy > maxEt )
         maxEt = et;

      // Stop phi wrap
      float dPhi1 = etaphiCorners[2].fY - etaphiCorners[1].fY;
      float dPhi2 = etaphiCorners[3].fY - etaphiCorners[0].fY;
      float dPhi3 = etaphiCorners[1].fY - etaphiCorners[2].fY;
      float dPhi4 = etaphiCorners[0].fY - etaphiCorners[3].fY;

      if( dPhi1 > 1 )
         etaphiCorners[2].fY = etaphiCorners[2].fY - ( 2 * TMath::Pi() );
      if( dPhi2 > 1 )
         etaphiCorners[3].fY = etaphiCorners[3].fY - ( 2 * TMath::Pi() );
      if( dPhi3 > 1 )
         etaphiCorners[2].fY = etaphiCorners[2].fY + ( 2 * TMath::Pi() );
      if( dPhi4 > 1 )
         etaphiCorners[3].fY = etaphiCorners[3].fY + ( 2 * TMath::Pi() );

      FWPFLegoRecHit *recHit = new FWPFLegoRecHit( etaphiCorners, itemHolder, this, vc, energy, et );
      recHit->setSquareColor(item()->defaultDisplayProperties().color());
      m_recHits.push_back( recHit );
   }
      
      m_maxEnergy = maxEnergy;
      m_maxEt = maxEt;
      m_maxEnergyLog = log(maxEnergy);
      m_maxEtLog = log(maxEt);

      scaleProduct( product, FWViewType::kLegoPFECAL, vc );
}

//______________________________________________________________________________________________________
void
FWPFEcalRecHitLegoProxyBuilder::cleanLocal()
{
   for( std::vector<FWPFLegoRecHit*>::iterator i = m_recHits.begin(); i != m_recHits.end(); ++i )
      delete (*i);

   m_recHits.clear();
}

//______________________________________________________________________________________________________
REGISTER_FWPROXYBUILDER( FWPFEcalRecHitLegoProxyBuilder, EcalRecHit, "Ecal RecHit", FWViewType::kLegoPFECALBit );
