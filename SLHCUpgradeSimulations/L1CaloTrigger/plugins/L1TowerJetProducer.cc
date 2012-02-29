
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/L1CaloAlgoBase.h"

#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "SimDataFormats/SLHC/interface/L1CaloTowerFwd.h"

#include <algorithm> 
#include <string> 
#include <vector>
 

class L1TowerJetProducer:
public L1CaloAlgoBase < l1slhc::L1CaloTowerCollection, l1slhc::L1TowerJetCollection >
{
  public:
	L1TowerJetProducer( const edm::ParameterSet & );
	 ~L1TowerJetProducer(  );

	// void initialize( );

	void algorithm( const int &, const int & );

  private:
	void calculateJetPosition( l1slhc::L1TowerJet & lJet );

	//some helpful members
	int mJetSize;
	l1slhc::L1TowerJet::tJetShape mJetShape;

	std::vector< std::pair< int , int > > mJetShapeMap;

};

L1TowerJetProducer::L1TowerJetProducer( const edm::ParameterSet & aConfig ):L1CaloAlgoBase < l1slhc::L1CaloTowerCollection, l1slhc::L1TowerJetCollection > ( aConfig )
{
	mJetSize = aConfig.getParameter<unsigned>("JetSize");
	mPhiOffset = -mJetSize;
	mEtaOffset = -mJetSize;
	// mPhiIncrement = 1;
	// mEtaIncrement = 1;

	mJetShapeMap.reserve(256); //jets will never be 16x16 but it is a nice round number

	std::string lJetShape = aConfig.getParameter< std::string >("JetShape");
	std::transform( lJetShape.begin() , lJetShape.end() , lJetShape.begin() , ::toupper ); //do the comparison in upper case so config file can read "Circle", "circle", "CIRCLE", "cIrClE", etc. and give the same result.

	std::cout << "Creating JetShapeMap:" << std::endl;

	if ( lJetShape == "CIRCLE" ){
		mJetShape = l1slhc::L1TowerJet::circle;


		double lCentre( (mJetSize-1) / 2.0 );
		double lDelta;

		std::vector<double> lDeltaSquare;
		for( int i = 0 ; i != mJetSize ; ++i ){
			lDelta = double(i) - lCentre;
			lDeltaSquare.push_back( lDelta*lDelta );
		}

		double lDeltaRSquare;
		double lDeltaRSquareMax( (mJetSize*mJetSize) / 4.0 );

		for( int x = 0 ; x != mJetSize ; ++x ){
			for( int y = 0 ; y != mJetSize ; ++y ){
				lDeltaRSquare = lDeltaSquare[x] + lDeltaSquare[y];
				if( lDeltaRSquare <= lDeltaRSquareMax ){
					mJetShapeMap.push_back( std::make_pair( x , y ) );
					std::cout << "#" << std::flush;
				}else{
					std::cout << "-" << std::flush;
				}
			}
			std::cout << std::endl;
		}

	}else{

		mJetShape = l1slhc::L1TowerJet::square;

		for( int x = 0 ; x != mJetSize ; ++x ){
			for( int y = 0 ; y != mJetSize ; ++y ){
				mJetShapeMap.push_back( std::make_pair( x , y ) );
			}
		}
	}

	std::cout << "JetShapeMap includes " << mJetShapeMap.size() << " towers." << std::endl;

}

L1TowerJetProducer::~L1TowerJetProducer(  )
{
}

/* 
   void L1TowerJetProducer::initialize( ) { }
*/

void L1TowerJetProducer::algorithm( const int &aEta, const int &aPhi )
{


	int lTowerIndex = mCaloTriggerSetup->getBin( aEta, aPhi );
	std::pair < int, int > lTowerEtaPhi = mCaloTriggerSetup->getTowerEtaPhi( lTowerIndex );


	l1slhc::L1TowerJet lJet( mJetSize, mJetShape , lTowerEtaPhi.first , lTowerEtaPhi.second );

	for ( std::vector< std::pair< int , int > >::const_iterator lJetShapeMapIt = mJetShapeMap.begin() ; lJetShapeMapIt != mJetShapeMap.end() ; ++lJetShapeMapIt )
	{
		l1slhc::L1CaloTowerCollection::const_iterator lTowerItr = fetch( aEta+(lJetShapeMapIt->first) , aPhi+(lJetShapeMapIt->second) );
		if ( lTowerItr != mInputCollection->end(  ) )
		{
			l1slhc::L1CaloTowerRef lRef( mInputCollection, lTowerItr - mInputCollection->begin(  ) );
			lJet.addConstituent( lRef );
		}
	}

	if ( lJet.E(  ) > 0 )
	{
		calculateJetPosition( lJet );
		mOutputCollection->insert( lTowerEtaPhi.first, lTowerEtaPhi.second, lJet );
		//std::cout << "towerjet " << lJet.p4(  ) << std::endl;
	}

}





// This function needs validating!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void L1TowerJetProducer::calculateJetPosition( l1slhc::L1TowerJet & lJet )
{

	double eta;		
	double halfTowerOffset = 0.0435;

	double halfJetSize = double(lJet.JetSize()) / 2.0;

	int abs_eta = abs( lJet.iEta(  ) + int(halfJetSize) );

	if ( abs_eta < 21 )
	{
//		eta = ( abs_eta * 0.0870 ) - halfTowerOffset;
//		if( lJet.JetSize() % 2 == 1 ){
//			eta += halfTowerOffset;
//		}

		eta = ( abs_eta * 0.0870 );
		if( lJet.JetSize() % 2 == 0 ){
			eta -= halfTowerOffset;
		}

	}
	else
	{
		const double endcapEta[8] = { 0.09, 0.1, 0.113, 0.129, 0.15, 0.178, 0.15, 0.35 };
		abs_eta -= 21;
	
		eta = 1.74;

//		for ( int i = 0; i <= abs_eta; ++i )
//		{
//			eta += endcapEta[i];
//		}
//		eta -= endcapEta[abs_eta] / 2.;

		for ( int i = 0; i != abs_eta; ++i )
		{
			eta += endcapEta[i];
		}
		if( lJet.JetSize() % 2 == 0 ){
			eta += endcapEta[abs_eta] / 2.;
		}


	}

	if ( lJet.iEta(  ) < -halfJetSize )
		eta = -eta;

	double phi = ( ( lJet.iPhi(  ) + halfJetSize ) * 0.087 ) - halfTowerOffset;
	double Et = double( lJet.E(  ) ) / 2.;

	lJet.setP4( math::PtEtaPhiMLorentzVector( Et, eta, phi, 0. ) );
}



DEFINE_EDM_PLUGIN( edm::MakerPluginFactory, edm::WorkerMaker < L1TowerJetProducer >, "L1TowerJetProducer" );
DEFINE_FWK_PSET_DESC_FILLER( L1TowerJetProducer );
