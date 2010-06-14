#include "TEveWindow.h"

#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/Core/interface/FWColorManager.h"
#include "Fireworks/Core/interface/FWL1TriggerTableView.h"
#include "Fireworks/Core/interface/FWL1TriggerTableViewManager.h"
#include "Fireworks/Core/interface/FWL1TriggerTableViewTableManager.h"
#include "Fireworks/Core/interface/fwLog.h"
#include "Fireworks/TableWidget/interface/FWTableWidget.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtTriggerMenuLite.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <iomanip>

static const std::string kTableView = "L1TriggerTableView";
static const std::string kColumns = "columns";
static const std::string kSortColumn = "sortColumn";
static const std::string kDescendingSort = "descendingSort";

FWL1TriggerTableView::FWL1TriggerTableView( TEveWindowSlot* parent, FWL1TriggerTableViewManager *manager )
   : m_manager( manager ),
     m_tableManager( new FWL1TriggerTableViewTableManager( this ) ),
     m_tableWidget( 0 ),
     m_currentColumn( -1 )
{
	m_columns.push_back( Column( "Algorithm Name" ));
	m_columns.push_back( Column( "Result" ) );
	m_columns.push_back( Column( "Bit Number" ) );
	m_columns.push_back( Column( "Prescale" ) );
   m_eveWindow = parent->MakeFrame( 0 );
   TGCompositeFrame *frame = m_eveWindow->GetGUICompositeFrame();

   m_vert = new TGVerticalFrame( frame );
   frame->AddFrame( m_vert, new TGLayoutHints( kLHintsExpandX | kLHintsExpandY ));

   m_tableWidget = new FWTableWidget( m_tableManager, m_vert );
   resetColors( m_manager->colorManager() );
   m_tableWidget->SetHeaderBackgroundColor( gVirtualX->GetPixel( kWhite ));
   m_tableWidget->Connect( "columnClicked(Int_t,Int_t,Int_t)", "FWL1TriggerTableView",
								  this, "columnSelected(Int_t,Int_t,Int_t)");
   m_vert->AddFrame( m_tableWidget, new TGLayoutHints( kLHintsExpandX | kLHintsExpandY ));
   dataChanged();
   frame->MapSubwindows();
   frame->Layout();
   frame->MapWindow();
}

FWL1TriggerTableView::~FWL1TriggerTableView( void )
{
   // take out composite frame and delete it directly (without the timeout)
   TGCompositeFrame *frame = m_eveWindow->GetGUICompositeFrame();
   frame->RemoveFrame( m_vert );
   delete m_vert;

   m_eveWindow->DestroyWindowAndSlot();
   delete m_tableManager;
}

void
FWL1TriggerTableView::setBackgroundColor( Color_t iColor )
{
	m_tableWidget->SetBackgroundColor( gVirtualX->GetPixel( iColor ));
}

void FWL1TriggerTableView::resetColors( const FWColorManager &manager )
{
   m_tableWidget->SetBackgroundColor( gVirtualX->GetPixel( manager.background()));
   m_tableWidget->SetLineSeparatorColor( gVirtualX->GetPixel( manager.foreground()));
}

TGFrame*
FWL1TriggerTableView::frame( void ) const
{
   return 0;
}

const std::string&
FWL1TriggerTableView::typeName( void ) const
{
   return staticTypeName();
}

void
FWL1TriggerTableView::saveImageTo( const std::string& iName ) const
{
	fwLog( fwlog::kWarning ) << "FWL1TriggerTableView::saveImageTo is not implemented." << std::endl;
}

void FWL1TriggerTableView::dataChanged( void )
{
   m_columns.at(0).values.clear();
   m_columns.at(1).values.clear();
   m_columns.at(2).values.clear();
   m_columns.at(3).values.clear();
   if(! m_manager->items().empty() && m_manager->items().front() != 0 )
   {
		if( fwlite::Event* event = const_cast<fwlite::Event*>( m_manager->items().front()->getEvent()))
      {
			fwlite::Handle<L1GtTriggerMenuLite> triggerMenuLite;
			fwlite::Handle<L1GlobalTriggerReadoutRecord> triggerRecord;

			try
			{
				// FIXME: Replace magic strings with configurable ones
				triggerMenuLite.getByLabel( event->getRun(), "l1GtTriggerMenuLite", "", "" );
				triggerRecord.getByLabel( *event, "gtDigis", "", "" );
			}
			catch( cms::Exception& )
			{
				fwLog( fwlog::kWarning ) << "FWL1TriggerTableView: no L1Trigger menu is available." << std::endl;
				m_tableManager->dataChanged();
				return;
			}
	  
			if( triggerMenuLite.isValid() && triggerRecord.isValid() )
			{
				const L1GtTriggerMenuLite::L1TriggerMap& algorithmMap = triggerMenuLite->gtAlgorithmMap();
				
				int pfIndexTechTrig = -1;
				int pfIndexAlgoTrig = -1;

				/// prescale factors
				std::vector<std::vector<int> > prescaleFactorsAlgoTrig = triggerMenuLite->gtPrescaleFactorsAlgoTrig();
				std::vector<std::vector<int> > prescaleFactorsTechTrig = triggerMenuLite->gtPrescaleFactorsTechTrig();
				pfIndexAlgoTrig = ( triggerRecord->gtFdlWord()).gtPrescaleFactorIndexAlgo();
				pfIndexTechTrig = ( triggerRecord->gtFdlWord()).gtPrescaleFactorIndexTech();

				const DecisionWord dWord = triggerRecord->decisionWord();
				for( L1GtTriggerMenuLite::CItL1Trig itTrig = algorithmMap.begin(), itTrigEnd = algorithmMap.end();
					 itTrig != itTrigEnd; ++itTrig )
				{
					const unsigned int bitNumber = itTrig->first;
					const std::string& aName = itTrig->second;
					int errorCode = 0;
					const bool result = triggerMenuLite->gtTriggerResult( aName, dWord, errorCode );

					m_columns.at(0).values.push_back( aName );
					m_columns.at(1).values.push_back( Form( "%d", result ));
					m_columns.at(2).values.push_back( Form( "%d", bitNumber ));
					m_columns.at(3).values.push_back( Form( "%d", prescaleFactorsAlgoTrig.at( pfIndexAlgoTrig ).at( bitNumber )));
				}
				const TechnicalTriggerWord ttWord = triggerRecord->technicalTriggerWord();
				
				int tBitNumber = 0;
				int tBitResult = 0;
				for( TechnicalTriggerWord::const_iterator tBitIt = ttWord.begin(), tBitEnd = ttWord.end(); 
					 tBitIt != tBitEnd; ++tBitIt, ++tBitNumber )
				{
					if( *tBitIt )
						tBitResult = 1;
					else
						tBitResult = 0;

					m_columns.at(0).values.push_back( "TechTrigger" );
					m_columns.at(1).values.push_back( Form( "%d", tBitResult ));
					m_columns.at(2).values.push_back( Form( "%d", tBitNumber ));
					m_columns.at(3).values.push_back( Form( "%d", prescaleFactorsTechTrig.at( pfIndexTechTrig ).at( tBitNumber )));
				}
			}
         else
         {
				fwLog( fwlog::kWarning ) << "FWL1TriggerTableView: " <<
				"L1GtTriggerMenuLite.isValid()=" << triggerMenuLite.isValid() << ", " <<
				"L1GlobalTriggerReadoutRecord.isValid()=" << triggerRecord.isValid() << "." << std::endl;
         }
      }
   }
   
   m_tableManager->dataChanged();
}

void
FWL1TriggerTableView::columnSelected( Int_t iCol, Int_t iButton, Int_t iKeyMod )
{
   if( iButton == 1 || iButton == 3 )
      m_currentColumn = iCol;
}

// void 
// FWL1TriggerTableView::updateFilter(void)
// {
//    dataChanged();
// }

//
// static member functions
//
const std::string&
FWL1TriggerTableView::staticTypeName( void )
{
   static std::string s_name( "L1TriggerTable" );
   return s_name;
}

