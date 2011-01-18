#include "HLTriggerOffline/Egamma/interface/EmDQMPostProcessor.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include <iostream>
#include <string.h>
#include <iomanip>
#include <fstream>
#include <math.h>


//----------------------------------------------------------------------

EmDQMPostProcessor::EmDQMPostProcessor(const edm::ParameterSet& pset)
{
  subDir_  = pset.getUntrackedParameter<std::string>("subDir");

  dataSet_ = pset.getUntrackedParameter<std::string>("dataSet","unknown");

  normalizeToReco = pset.getUntrackedParameter<bool>("normalizeToReco",false);
}

//----------------------------------------------------------------------

void EmDQMPostProcessor::endRun(edm::Run const& run, edm::EventSetup const& es)
{
  //////////////////////////////////
  // setup DQM stor               //
  //////////////////////////////////
  
  DQMStore * dqm = 0;
  dqm = edm::Service<DQMStore>().operator->();

  if ( ! dqm ) {
    edm::LogInfo("EmDQMPostProcessor") << "Cannot create DQMStore instance\n";
    return;
  }


  //go to the directory to be processed
  if(dqm->dirExists(subDir_)) dqm->cd(subDir_);
  else {
    edm::LogWarning("EmDQMPostProcessor") << "cannot find directory: " << subDir_ << " , skipping";
    return;
  }

  //--------------------
  // with respect to what are (some) efficiencies calculated ?
  std::string shortReferenceName;
  if (normalizeToReco)
    shortReferenceName = "RECO";
  else
    shortReferenceName = "gen";



  //--------------------


  //////////////////////////////////
  //loop over all triggers/samples//
  //////////////////////////////////

  // store dataset name (if defined) in output file
  // DQMStore:bookString seems to put a key in the file which is
  // of the form <dataSet>s=....</dataSet> which is not very convenient
  // (it points to a null pointer, one would have to loop over
  // all keys of the corresponding directory in the ROOT file
  // and check whether it is of the desired form and then parse
  // it from this string...).
  //
  // So we store the name of the dataset as the title of a histogram,
  // which is much easier to access...
  // TH1D *dataSetNameHisto = 
  dqm->book1D("DataSetNameHistogram",dataSet_,1,0,1);

  std::vector<std::string> subdirectories = dqm->getSubdirs();
  for(std::vector<std::string>::iterator dir = subdirectories.begin() ;dir!= subdirectories.end(); dir++ ){
    dqm->cd(*dir);


    ////////////////////////////////////////////////////////
    // Do everything twice: once for mc-matched histos,   //
    // once for unmatched histos                          //
    ////////////////////////////////////////////////////////

    std::vector<std::string> postfixes;
    postfixes.push_back("");   //unmatched histograms
    postfixes.push_back("_RECO_matched"); // for data

    // we put this on the list even when we're running on 
    // data (where there is no generator information).
    // The first test in the loop will then fail and
    // the iteration is skipped.
    postfixes.push_back("_MC_matched"); 

    for(std::vector<std::string>::iterator postfix=postfixes.begin(); postfix!=postfixes.end();postfix++){
      
      /////////////////////////////////////
      // computer per-event efficiencies //
      /////////////////////////////////////
      
      std::string histoName="efficiency_by_step"+ *postfix;
      std::string baseName = "total_eff"+ *postfix;

      TH1F* basehist = getHistogram(dqm, dqm->pwd() + "/" + baseName);
      if (basehist == NULL)
	{
	  edm::LogWarning("EmDQMPostProcessor") << "histogram " << (dqm->pwd() + "/" + baseName) << " does not exist, skipping postfix '" << *postfix << "'"; 
	  continue;
	}
      TProfile* total = dqm->bookProfile(histoName,histoName,basehist->GetXaxis()->GetNbins(),basehist->GetXaxis()->GetXmin(),basehist->GetXaxis()->GetXmax(),0.,1.2)->getTProfile();
      total->GetXaxis()->SetBinLabel(1,basehist->GetXaxis()->GetBinLabel(1));
      
//       std::vector<std::string> mes = dqm->getMEs();
//       for(std::vector<std::string>::iterator me = mes.begin() ;me!= mes.end(); me++ )
// 	std::cout <<*me <<std::endl;
//       std::cout <<std::endl;
      
      double value=0;
      double errorh=0,errorl=0,error=0;    
      //compute stepwise total efficiencies 
      for(int bin= total->GetNbinsX()-2 ; bin > 1  ; bin--){
	value=0;
	errorl=0;
	errorh=0;
	error=0;
	if(basehist->GetBinContent(bin-1) != 0){
	  Efficiency( (int)basehist->GetBinContent(bin), (int)basehist->GetBinContent(bin-1), 0.683, value, errorl, errorh );
	  error = value-errorl>errorh-value ? value-errorl : errorh-value;
	}
	total->SetBinContent( bin, value );
	total->SetBinEntries( bin, 1 );
	total->SetBinError( bin, sqrt(value*value+error*error) );
	total->GetXaxis()->SetBinLabel(bin,basehist->GetXaxis()->GetBinLabel(bin));
      }

      //set first bin to L1 efficiency
      if(basehist->GetBinContent(basehist->GetNbinsX()) !=0 ){
	Efficiency( (int)basehist->GetBinContent(1), (int)basehist->GetBinContent(basehist->GetNbinsX()), 0.683, value, errorl, errorh );
	error= value-errorl>errorh-value ? value-errorl : errorh-value;
      }else{
	value=0;error=0;
      }
      total->SetBinContent(1,value);
      total->SetBinEntries(1, 1 );
      total->SetBinError(1, sqrt(value*value+error*error) );
     
      // total efficiency relative to gen or reco
      if(basehist->GetBinContent(basehist->GetNbinsX()) !=0 ){
	Efficiency( (int)basehist->GetBinContent(basehist->GetNbinsX()-2), (int)basehist->GetBinContent(basehist->GetNbinsX()), 0.683, value, errorl, errorh );
	error= value-errorl>errorh-value ? value-errorl : errorh-value;
      }else{
	value=0;error=0;
      }
      total->SetBinContent(total->GetNbinsX(),value);
      total->SetBinEntries(total->GetNbinsX(),1);
      total->SetBinError(total->GetNbinsX(),sqrt(value*value+error*error));
      total->GetXaxis()->SetBinLabel(total->GetNbinsX(),("total efficiency rel. " + shortReferenceName).c_str());
      
      //total efficiency relative to L1
      if(basehist->GetBinContent(1) !=0 ){
	Efficiency( (int)basehist->GetBinContent(basehist->GetNbinsX()-2), (int)basehist->GetBinContent(1), 0.683, value, errorl, errorh );
	error= value-errorl > errorh-value ? value-errorl : errorh-value;
      }else{
	value=0;error=0;
      }
      total->SetBinContent(total->GetNbinsX()-1,value);
      total->SetBinError(total->GetNbinsX()-1,sqrt(value*value+error*error));
      total->SetBinEntries(total->GetNbinsX()-1,1);
      total->GetXaxis()->SetBinLabel(total->GetNbinsX()-1,"total efficiency rel. L1");

      //----------------------------------------

      ///////////////////////////////////////////
      // compute per-object efficiencies       //
      ///////////////////////////////////////////
      //MonitorElement *eff, *num, *denom, *genPlot, *effVsGen, *effL1VsGen;
      std::vector<std::string> varNames; 
      varNames.push_back("eta"); 
      varNames.push_back("phi"); 
      varNames.push_back("et");

      std::string filterName;
      std::string filterName2;
      std::string denomName;
      std::string numName;

      // Get the gen-level (or reco, for data) plots
      std::string genName;

      // Get the L1 over gen filter first
      filterName2= total->GetXaxis()->GetBinLabel(1);
	
      //loop over variables (eta/phi/et)
      for(std::vector<std::string>::iterator var = varNames.begin(); var != varNames.end() ; var++){
	
	numName   = dqm->pwd() + "/" + filterName2 + *var + *postfix;

	if (normalizeToReco)
	  genName   = dqm->pwd() + "/reco_" + *var ;
	else
	  genName   = dqm->pwd() + "/gen_" + *var ;

	// Create the efficiency plot
	if(!dividehistos(dqm,numName,genName,"efficiency_"+filterName2+"_vs_"+*var +*postfix,*var,"eff. of"+filterName2+" vs "+*var +*postfix))
	  break;
      } // loop over variables
    
      // get the filter names from the bin-labels of the master-histogram
      for (int filter=1; filter < total->GetNbinsX()-2; filter++) {
	filterName = total->GetXaxis()->GetBinLabel(filter);
	filterName2= total->GetXaxis()->GetBinLabel(filter+1);
	
	//loop over variables (eta/et/phi)
	for(std::vector<std::string>::iterator var = varNames.begin(); var != varNames.end() ; var++){
	  numName   = dqm->pwd() + "/" + filterName2 + *var + *postfix;
	  denomName = dqm->pwd() + "/" + filterName  + *var + *postfix;

	  // Is this the last filter? Book efficiency vs gen (or reco, for data) level
	  std::string temp = *postfix;
          if (filter==total->GetNbinsX()-3 && temp.find("matched")!=std::string::npos) {
	    if (normalizeToReco)
	      genName = dqm->pwd() + "/reco_" + *var;
	    else
	      genName = dqm->pwd() + "/gen_" + *var;

	    if(!dividehistos(dqm,numName,genName,"final_eff_vs_"+*var,*var,"Efficiency Compared to " + shortReferenceName + " vs "+*var))
	      break;
	  }

	  if(!dividehistos(dqm,numName,denomName,"efficiency_"+filterName2+"_vs_"+*var +*postfix,*var,"efficiency_"+filterName2+"_vs_"+*var + *postfix))
	    break;

	} // loop over variables
      } // loop over monitoring modules within path
    } // loop over postfixes 
    dqm->goUp();
  }
  
}

//----------------------------------------------------------------------

TProfile* EmDQMPostProcessor::dividehistos(DQMStore * dqm, const std::string& numName, const std::string& denomName, const std::string& outName,const std::string& label,const std::string& titel){
  //std::cout << numName <<std::endl;

  TH1F* num = getHistogram(dqm,numName);

  //std::cout << denomName << std::endl;
  TH1F* denom = getHistogram(dqm, denomName);

  if (num == NULL)
    edm::LogWarning("EmDQMPostProcessor") << "numerator histogram " << numName << " does not exist";

  if (denom == NULL)
    edm::LogWarning("EmDQMPostProcessor") << "denominator histogram " << denomName << " does not exist";

  // Check if histograms actually exist

  if(!num || !denom) return 0;

  // Make sure we are able to book new element
  if (!dqm) return 0;
  
  TProfile* out = dqm->bookProfile(outName,titel,num->GetXaxis()->GetNbins(),num->GetXaxis()->GetXmin(),num->GetXaxis()->GetXmax(),0.,1.2)->getTProfile();
  out->GetXaxis()->SetTitle(label.c_str());
  out->SetYTitle("Efficiency");
  out->SetOption("PE");
  out->SetLineColor(2);
  out->SetLineWidth(2);
  out->SetMarkerStyle(20);
  out->SetMarkerSize(0.8);
  out->SetStats(kFALSE);
  for(int i=1;i<=num->GetNbinsX();i++){
    double e, low, high;
    Efficiency( (int)num->GetBinContent(i), (int)denom->GetBinContent(i), 0.683, e, low, high );
    double err = e-low>high-e ? e-low : high-e;
    //here is the trick to store info in TProfile:
    out->SetBinContent( i, e );
    out->SetBinEntries( i, 1 );
    out->SetBinError( i, sqrt(e*e+err*err) );
  }

  return out;
}

//----------------------------------------------------------------------

TH1F *
EmDQMPostProcessor::getHistogram(DQMStore *dqm, const std::string &histoPath)
{
  MonitorElement *monElement = dqm->get(histoPath);
  if (monElement != NULL)
    return monElement->getTH1F();
  else
    return NULL;
}

//----------------------------------------------------------------------

DEFINE_FWK_MODULE(EmDQMPostProcessor);

//----------------------------------------------------------------------
