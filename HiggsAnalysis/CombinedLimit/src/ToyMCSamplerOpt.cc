#include "../interface/ToyMCSamplerOpt.h"
#include "../interface/utils.h"
#include <memory>
#include <stdexcept>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <RooSimultaneous.h>
#include <RooRealVar.h>
#include <RooProdPdf.h>
#include <RooPoisson.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooRandom.h>

ToyMCSamplerOpt::ToyMCSamplerOpt(RooStats::TestStatistic& ts, Int_t ntoys, RooAbsPdf *globalObsPdf) :
    ToyMCSampler(ts, ntoys),
    globalObsPdf_(globalObsPdf),
    globalObsValues_(0), globalObsIndex_(-1),
    weightVar_(0),
    _allVars(0)
{
}


ToyMCSamplerOpt::ToyMCSamplerOpt(const RooStats::ToyMCSampler &base) :
    ToyMCSampler(base),
    globalObsPdf_(0),
    globalObsValues_(0), globalObsIndex_(-1),
    weightVar_(0),
    _allVars(0)
{
}

ToyMCSamplerOpt::ToyMCSamplerOpt(const ToyMCSamplerOpt &other) :
    ToyMCSampler(other),
    globalObsPdf_(0),
    globalObsValues_(0), globalObsIndex_(-1),
    weightVar_(0),
    _allVars(0)
{
}

ToyMCSamplerOpt::~ToyMCSamplerOpt()
{
    delete weightVar_;
    for (std::map<RooAbsPdf *, toymcoptutils::SimPdfGenInfo *>::iterator it = genCache_.begin(), ed = genCache_.end(); it != ed; ++it) {
        delete it->second;
    }
    genCache_.clear();
    delete _allVars;
    delete globalObsValues_;
}


toymcoptutils::SinglePdfGenInfo::SinglePdfGenInfo(RooAbsPdf &pdf, const RooArgSet& observables, bool preferBinned, const RooDataSet* protoData, int forceEvents) :
   mode_(pdf.canBeExtended() ? (preferBinned ? Binned : Unbinned) : Counting),
   pdf_(&pdf),
   spec_(0),weightVar_(0)
{
   if (pdf.canBeExtended()) {
       if (pdf.getAttribute("forceGenBinned")) mode_ = Binned;
       else if (pdf.getAttribute("forceGenPoisson")) mode_ = Poisson;
       else if (pdf.getAttribute("forceGenUnbinned")) mode_ = Unbinned;
       //else std::cout << "Pdf " << pdf.GetName() << " has no preference" << std::endl;
   }

   RooArgSet *obs = pdf.getObservables(observables);
   observables_.add(*obs);
   delete obs;
   //if (mode_ == Unbinned) spec_ = protoData ? pdf.prepareMultiGen(observables_, RooFit::Extended(), RooFit::ProtoData(*protoData, true, true)) 
   //                                         : pdf.prepareMultiGen(observables_, RooFit::Extended());
}

toymcoptutils::SinglePdfGenInfo::~SinglePdfGenInfo()
{
    delete spec_;
    delete weightVar_;
}


RooAbsData *  
toymcoptutils::SinglePdfGenInfo::generate(const RooDataSet* protoData, int forceEvents) 
{
    assert(forceEvents == 0 && "SinglePdfGenInfo: forceEvents must be zero at least for now");
    RooAbsData *ret = 0;
    if (mode_ == Unbinned) {
        //ret = pdf_->generate(*spec_);
        ret = pdf_->generate(observables_, RooFit::Extended());
    } else if (mode_ == Binned) {
        //ret = protoData ? pdf_->generateBinned(observables_, RooFit::Extended(), RooFit::ProtoData(*protoData, true, true))
        //                : pdf_->generateBinned(observables_, RooFit::Extended());
        // generateBinnedWorkaround
        RooDataSet *data =  pdf_->generate(observables_, RooFit::Extended());
        ret = new RooDataHist(data->GetName(), "", *data->get(), *data);
        delete data;
    } else if (mode_ == Poisson) {
        return generateWithHisto(weightVar_, false);
    } else if (mode_ == Counting) {
        ret = pdf_->generate(observables_, 1);
    } else throw std::logic_error("Mode not foreseen in SinglePdfGenInfo::generate");
    //std::cout << "Dataset generated from " << pdf_->GetName() << " (weighted? " << ret->isWeighted() << ")" << std::endl;
    //utils::printRAD(ret);
    return ret;
}

RooDataSet *  
toymcoptutils::SinglePdfGenInfo::generateAsimov(RooRealVar *&weightVar) 
{
    return generateWithHisto(weightVar, true);
}

RooDataSet *  
toymcoptutils::SinglePdfGenInfo::generateWithHisto(RooRealVar *&weightVar, bool asimov) 
{
    if (mode_ == Counting) return generateCountingAsimov();
    if (observables_.getSize() > 3) throw std::invalid_argument(std::string("ERROR in SinglePdfGenInfo::generateWithHisto for ") + pdf_->GetName() + ", more than 3 observable");
    RooArgList obs(observables_);
    RooRealVar *x = (RooRealVar*)obs.at(0);
    RooRealVar *y = obs.getSize() > 1 ? (RooRealVar*)obs.at(1) : 0;
    RooRealVar *z = obs.getSize() > 2 ? (RooRealVar*)obs.at(2) : 0;
    if (weightVar == 0) weightVar = new RooRealVar("_weight_","",1.0);

    RooCmdArg ay = (y ? RooFit::YVar(*y) : RooCmdArg::none());
    RooCmdArg az = (z ? RooFit::YVar(*z) : RooCmdArg::none());
    std::auto_ptr<TH1> hist(pdf_->createHistogram("htemp", *x, ay, az));

    double expectedEvents = pdf_->expectedEvents(observables_);
    hist->Scale(expectedEvents/ hist->Integral()); 
    RooArgSet obsPlusW(obs); obsPlusW.add(*weightVar);
    RooDataSet *data = new RooDataSet(TString::Format("%sData", pdf_->GetName()), "", obsPlusW, weightVar->GetName());
    switch (obs.getSize()) {
        case 1:
            for (int i = 1, n = hist->GetNbinsX(); i <= n; ++i) {
                x->setVal(hist->GetXaxis()->GetBinCenter(i));
                data->add(observables_, asimov ? hist->GetBinContent(i) : RooRandom::randomGenerator()->Poisson(hist->GetBinContent(i)) );
            }
            break;
        case 2:
            {
            TH2& h2 = dynamic_cast<TH2&>(*hist);
            for (int ix = 1, nx = h2.GetNbinsX(); ix <= nx; ++ix) {
            for (int iy = 1, ny = h2.GetNbinsY(); iy <= ny; ++iy) {
                x->setVal(h2.GetXaxis()->GetBinCenter(ix));
                y->setVal(h2.GetYaxis()->GetBinCenter(iy));
                data->add(observables_, asimov ? h2.GetBinContent(ix,iy) : RooRandom::randomGenerator()->Poisson(h2.GetBinContent(ix,iy)) );
            } }
            }
            break;
        case 3:
            {
            TH3& h3 = dynamic_cast<TH3&>(*hist);
            for (int ix = 1, nx = h3.GetNbinsX(); ix <= nx; ++ix) {
            for (int iy = 1, ny = h3.GetNbinsY(); iy <= ny; ++iy) {
            for (int iz = 1, nz = h3.GetNbinsZ(); iz <= nz; ++iz) {
                x->setVal(h3.GetXaxis()->GetBinCenter(ix));
                y->setVal(h3.GetYaxis()->GetBinCenter(iy));
                z->setVal(h3.GetYaxis()->GetBinCenter(iz));
                data->add(observables_, asimov ? h3.GetBinContent(ix,iy,iz) : RooRandom::randomGenerator()->Poisson(h3.GetBinContent(ix,iy,iz)) );
            } } }
            }
    }
    //std::cout << "Asimov dataset generated from " << pdf_->GetName() << " (sumw? " << data->sumEntries() << ", expected events " << expectedEvents << ")" << std::endl;
    //utils::printRDH(data);
    return data;
}


RooDataSet *  
toymcoptutils::SinglePdfGenInfo::generateCountingAsimov() 
{
    RooArgSet obs(observables_);
    RooProdPdf *prod = dynamic_cast<RooProdPdf *>(pdf_);
    RooPoisson *pois = 0;
    if (prod != 0) {
        setToExpected(*prod, observables_);
    } else if ((pois = dynamic_cast<RooPoisson *>(pdf_)) != 0) {
        setToExpected(*pois, observables_);
    } else throw std::logic_error("A counting model pdf must be either a RooProdPdf or a RooPoisson");
    RooDataSet *ret = new RooDataSet(TString::Format("%sData", pdf_->GetName()), "", obs);
    ret->add(obs);
    return ret;
}

void
toymcoptutils::SinglePdfGenInfo::setToExpected(RooProdPdf &prod, RooArgSet &obs) 
{
    std::auto_ptr<TIterator> iter(prod.pdfList().createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
        if (!a->dependsOn(obs)) continue;
        RooPoisson *pois = 0;
        if ((pois = dynamic_cast<RooPoisson *>(a)) != 0) {
            setToExpected(*pois, obs);
        } else {
            RooProdPdf *subprod = dynamic_cast<RooProdPdf *>(a);
            if (subprod) setToExpected(*subprod, obs);
            else throw std::logic_error("Illegal term in counting model: depends on observables, but not Poisson or Product");
        }
    }
}

void
toymcoptutils::SinglePdfGenInfo::setToExpected(RooPoisson &pois, RooArgSet &obs) 
{
    RooRealVar *myobs = 0;
    RooAbsReal *myexp = 0;
    std::auto_ptr<TIterator> iter(pois.serverIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
        if (obs.contains(*a)) {
            assert(myobs == 0 && "SinglePdfGenInfo::setToExpected(RooPoisson): Two observables??");
            myobs = dynamic_cast<RooRealVar *>(a);
            assert(myobs != 0 && "SinglePdfGenInfo::setToExpected(RooPoisson): Observables is not a RooRealVar??");
        } else {
            assert(myexp == 0 && "SinglePdfGenInfo::setToExpected(RooPoisson): Two expecteds??");
            myexp = dynamic_cast<RooAbsReal *>(a);
            assert(myexp != 0 && "SinglePdfGenInfo::setToExpected(RooPoisson): Expectedis not a RooAbsReal??");
        }
    }
    assert(myobs != 0 && "SinglePdfGenInfo::setToExpected(RooPoisson): No observable?");
    assert(myexp != 0 && "SinglePdfGenInfo::setToExpected(RooPoisson): No expected?");
    myobs->setVal(myexp->getVal());
}

toymcoptutils::SimPdfGenInfo::SimPdfGenInfo(RooAbsPdf &pdf, const RooArgSet& observables, bool preferBinned, const RooDataSet* protoData, int forceEvents) :
    pdf_(&pdf),
    cat_(0),
    observables_(observables),
    copyData_(true)
{
    assert(forceEvents == 0 && "SimPdfGenInfo: forceEvents must be zero at least for now");
    RooSimultaneous *simPdf = dynamic_cast<RooSimultaneous *>(&pdf);
    if (simPdf) {
        cat_ = const_cast<RooAbsCategoryLValue *>(&simPdf->indexCat());
        int nbins = cat_->numBins((const char *)0);
        pdfs_.resize(nbins, 0);
        RooArgList dummy;
        for (int ic = 0; ic < nbins; ++ic) {
            cat_->setBin(ic);
            RooAbsPdf *pdfi = simPdf->getPdf(cat_->getLabel());
            RooAbsPdf *newpdf = utils::factorizePdf(observables, *pdfi, dummy);
            pdfs_[ic] = new SinglePdfGenInfo(*newpdf, observables, preferBinned);
            if (newpdf != 0 && newpdf != pdfi) {
                ownedCrap_.addOwned(*newpdf); 
            }
        }
    } else {
        pdfs_.push_back(new SinglePdfGenInfo(pdf, observables, preferBinned, protoData, forceEvents));
    }
}

toymcoptutils::SimPdfGenInfo::~SimPdfGenInfo()
{
    for (std::vector<SinglePdfGenInfo *>::iterator it = pdfs_.begin(), ed = pdfs_.end(); it != ed; ++it) {
        delete *it;
    }
    pdfs_.clear();
    //for (std::map<std::string,RooDataSet*>::iterator it = datasetPieces_.begin(), ed = datasetPieces_.end(); it != ed; ++it) {
    for (std::map<std::string,RooAbsData*>::iterator it = datasetPieces_.begin(), ed = datasetPieces_.end(); it != ed; ++it) {
        delete it->second;
    }
    datasetPieces_.clear();
}


RooAbsData *  
toymcoptutils::SimPdfGenInfo::generate(RooRealVar *&weightVar, const RooDataSet* protoData, int forceEvents) 
{
    RooAbsData *ret = 0;
    TString retName =  TString::Format("%sData", pdf_->GetName());
    if (cat_ != 0) {
        //bool needsWeights = false;
        for (int i = 0, n = cat_->numBins((const char *)0); i < n; ++i) {
            if (pdfs_[i] == 0) continue;
            cat_->setBin(i);
            RooAbsData *&data =  datasetPieces_[cat_->getLabel()]; delete data;
            assert(protoData == 0);
            data = pdfs_[i]->generate(protoData); // I don't really know if protoData != 0 would make sense here
            if (data->isWeighted()) {
                if (weightVar == 0) weightVar = new RooRealVar("_weight_","",1.0);
                RooArgSet obs(*data->get()); 
                obs.add(*weightVar);
                RooDataSet *wdata = new RooDataSet(data->GetName(), "", obs, "_weight_");
                for (int i = 0, n = data->numEntries(); i < n; ++i) {
                    obs = *data->get(i);
                    if (data->weight()) wdata->add(obs, data->weight());
                }
                //std::cout << "DataHist was " << std::endl; utils::printRAD(data);
                delete data;
                data = wdata;
                //std::cout << "DataSet is " << std::endl; utils::printRAD(data);
            } 
            //if (data->isWeighted()) needsWeights = true;
        }
        if (copyData_) {
        // copyData is the "slow" mode used when generating toys with option "-t", 
        // that produces toys which can be saved standalone
#if 0 /// ===== Unfortunately this sometimes breaks the weights =======
            std::map<std::string,RooDataSet*> otherMap;
            for (std::map<std::string,RooAbsData*>::iterator it = datasetPieces_.begin(), ed = datasetPieces_.end(); it != ed; ++it) {
                RooDataSet* rds = dynamic_cast<RooDataSet*>(it->second);
                if (rds == 0) throw std::logic_error("Error, it should have been a RooDataSet");
                otherMap[it->first] = rds;
            }
            if (weightVar) {
                //std::cout << "Creating with weight" << std::endl;
                RooArgSet varsPlusWeight(observables_); varsPlusWeight.add(*weightVar);
                ret = new RooDataSet(retName, "", varsPlusWeight, RooFit::Index((RooCategory&)*cat_), RooFit::Import(otherMap), RooFit::WeightVar(*weightVar));
            } else {
                //std::cout << "Creating without weight" << std::endl;
                ret = new RooDataSet(retName, "", observables_, RooFit::Index((RooCategory&)*cat_), RooFit::Import(otherMap));
            }
#else //// ==== slower but safer solution
            RooArgSet vars(observables_), varsPlusWeight(observables_); 
            if (weightVar) varsPlusWeight.add(*weightVar);
            ret = new RooDataSet(retName, "", varsPlusWeight, (weightVar ? weightVar->GetName() : 0));
            for (std::map<std::string,RooAbsData*>::iterator it = datasetPieces_.begin(), ed = datasetPieces_.end(); it != ed; ++it) {
                cat_->setLabel(it->first.c_str());
                for (unsigned int i = 0, n = it->second->numEntries(); i < n; ++i) {
                    vars = *it->second->get(i);
                    ret->add(vars, it->second->weight());
                }
            }
#endif
        } else {
            // not copyData is the "fast" mode used when generating toys as a ToyMCSampler.
            // this doesn't copy the data, so the toys cannot outlive this class and each new
            // toy over-writes the memory of the previous one.
            ret = new RooDataSet(retName, "", observables_, RooFit::Index((RooCategory&)*cat_), RooFit::Link(datasetPieces_) /*, RooFit::OwnLinked()*/);
        }
    } else ret = pdfs_[0]->generate(protoData, forceEvents);
    //std::cout << "Dataset generated from sim pdf (weighted? " << ret->isWeighted() << ")" << std::endl; utils::printRAD(ret);
    return ret;
}

RooAbsData *  
toymcoptutils::SimPdfGenInfo::generateAsimov(RooRealVar *&weightVar) 
{
    RooAbsData *ret = 0;
    TString retName =  TString::Format("%sData", pdf_->GetName());
    if (cat_ != 0) {
        //bool needsWeights = false;
        for (int i = 0, n = cat_->numBins((const char *)0); i < n; ++i) {
            if (pdfs_[i] == 0) continue;
            cat_->setBin(i);
            RooAbsData *&data =  datasetPieces_[cat_->getLabel()]; delete data;
            data = pdfs_[i]->generateAsimov(weightVar); 
        }
        if (copyData_) { 
            // copyData is the "slow" mode used when generating toys with option "-t", 
            // that produces toys which can be saved standalone
            std::map<std::string,RooDataSet*> otherMap;
            for (std::map<std::string,RooAbsData*>::iterator it = datasetPieces_.begin(), ed = datasetPieces_.end(); it != ed; ++it) {
                RooDataSet* rds = dynamic_cast<RooDataSet*>(it->second);
                if (rds == 0) throw std::logic_error("Error, it should have been a RooDataSet");
                otherMap[it->first] = rds;
            }
            RooArgSet varsPlusWeight(observables_); varsPlusWeight.add(*weightVar);
            // here we can use the RooFit::Import because Asimov datasets are alwasy weighted 
            ret = new RooDataSet(retName, "", varsPlusWeight, RooFit::Index((RooCategory&)*cat_), RooFit::Import(otherMap), RooFit::WeightVar(*weightVar));
        } else {
            // not copyData is the "fast" mode used when generating toys as a ToyMCSampler.
            // this doesn't copy the data, so the toys cannot outlive this class and each new
            // toy over-writes the memory of the previous one.
            ret = new RooDataSet(retName, "", observables_, RooFit::Index((RooCategory&)*cat_), RooFit::Link(datasetPieces_) /*, RooFit::OwnLinked()*/);
        }
    } else ret = pdfs_[0]->generateAsimov(weightVar);
    //std::cout << "Asimov dataset generated from sim pdf " << pdf_->GetName() << " (sumw? " << ret->sumEntries() << ")" << std::endl; 
    //utils::printRAD(ret);
    return ret;
}


void
ToyMCSamplerOpt::SetPdf(RooAbsPdf& pdf) 
{
    ToyMCSampler::SetPdf(pdf);
    delete _allVars; _allVars = 0; 
    delete globalObsValues_; globalObsValues_ = 0; globalObsIndex_ = -1;
}

#if ROOT_VERSION_CODE < ROOT_VERSION(5,29,0)
//--- Taken from SVN HEAD ------------
RooAbsData* ToyMCSamplerOpt::GenerateToyData(RooArgSet& /*nullPOI*/) const {
   // This method generates a toy data set for the given parameter point taking
   // global observables into account.

   if (fObservables == NULL) { 
      ooccoutE((TObject*)NULL,InputArguments) << "Observables not set." << endl; 
      return 0; 
   }

   RooArgSet observables(*fObservables);
   if(fGlobalObservables  &&  fGlobalObservables->getSize()) {
      observables.remove(*fGlobalObservables);

      
      // generate one set of global observables and assign it
      // has problem for sim pdfs
      RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(fPdf);
      if(globalObsPdf_ || !simPdf){
        if (globalObsValues_ == 0 || globalObsIndex_ == globalObsValues_->numEntries()) {
            delete globalObsValues_;
            globalObsValues_ = (globalObsPdf_ ? globalObsPdf_ : fPdf)->generate(*fGlobalObservables, fNToys);
            globalObsIndex_  = 0;
        }
	const RooArgSet *values = globalObsValues_->get(globalObsIndex_++);
        //std::cout << "Generated for " << fPdf->GetName() << std::endl; values->Print("V");
	if (!_allVars) {
	  _allVars = fPdf->getObservables(*fGlobalObservables);
	}
	*_allVars = *values;

      } else {
#if 0
	if (_pdfList.size()==0) {
	  TIterator* citer = simPdf->indexCat().typeIterator() ;
	  RooCatType* tt = NULL;
	  while((tt=(RooCatType*) citer->Next())) {
	    RooAbsPdf* pdftmp = simPdf->getPdf(tt->GetName()) ;
	    RooArgSet* globtmp = pdftmp->getObservables(*fGlobalObservables) ;
	    RooAbsPdf::GenSpec* gs = pdftmp->prepareMultiGen(*globtmp,RooFit::NumEvents(1)) ;
	    _pdfList.push_back(pdftmp) ;
	    _obsList.push_back(globtmp) ;
	    _gsList.push_back(gs) ;
	  }
	}

	list<RooArgSet*>::iterator oiter = _obsList.begin() ;
	list<RooAbsPdf::GenSpec*>::iterator giter = _gsList.begin() ;
	for (list<RooAbsPdf*>::iterator iter = _pdfList.begin() ; iter != _pdfList.end() ; ++iter, ++giter, ++oiter) {
	  //RooDataSet* tmp = (*iter)->generate(**oiter,1) ;	  
	  RooDataSet* tmp = (*iter)->generate(**giter) ;
	  **oiter = *tmp->get(0) ;
	  delete tmp ;	  
	}	
#else
        //try fix for sim pdf
        TIterator* iter = simPdf->indexCat().typeIterator() ;
        RooCatType* tt = NULL;
        while((tt=(RooCatType*) iter->Next())) {

            // Get pdf associated with state from simpdf
            RooAbsPdf* pdftmp = simPdf->getPdf(tt->GetName()) ;

            // Generate only global variables defined by the pdf associated with this state
            RooArgSet* globtmp = pdftmp->getObservables(*fGlobalObservables) ;
            RooDataSet* tmp = pdftmp->generate(*globtmp,1) ;

            // Transfer values to output placeholder
            *globtmp = *tmp->get(0) ;

            // Cleanup
            delete globtmp ;
            delete tmp ;
        }
#endif
      } 
}

   RooAbsData* data = NULL;

   if(!fImportanceDensity) {
      // no Importance Sampling
      data = Generate(*fPdf, observables);
   }else{

      // Importance Sampling
      RooArgSet* allVars = fPdf->getVariables();
      RooArgSet* allVars2 = fImportanceDensity->getVariables();
      allVars->add(*allVars2);
      const RooArgSet* saveVars = (const RooArgSet*)allVars->snapshot();

      // the number of events generated is either the given fNEvents or
      // in case this is not given, the expected number of events of
      // the pdf with a Poisson fluctuation
      int forceEvents = 0;
      if(fNEvents == 0) {
         forceEvents = (int)fPdf->expectedEvents(observables);
         forceEvents = RooRandom::randomGenerator()->Poisson(forceEvents);
      }

      // need to be careful here not to overwrite the current state of the
      // nuisance parameters, ie they must not be part of the snapshot
      if(fImportanceSnapshot) *allVars = *fImportanceSnapshot;

      // generate with the parameters configured in this class
      //   NULL => no protoData
      //   overwriteEvents => replaces fNEvents it would usually take
      data = Generate(*fImportanceDensity, observables, NULL, forceEvents);

      *allVars = *saveVars;
      delete allVars;
      delete allVars2;
      delete saveVars;
   }

   return data;
}
#endif


RooAbsData *  
ToyMCSamplerOpt::Generate(RooAbsPdf& pdf, RooArgSet& observables, const RooDataSet* protoData, int forceEvents) const 
{
   if(fProtoData) {
      protoData = fProtoData;
      forceEvents = protoData->numEntries();
   }
   int events = forceEvents;
   if (events == 0) events = fNEvents;
   if (events != 0) return RooStats::ToyMCSampler::Generate(pdf, observables, protoData, forceEvents);
   toymcoptutils::SimPdfGenInfo *& info = genCache_[&pdf];
   if (info == 0) { 
       info = new toymcoptutils::SimPdfGenInfo(pdf, observables, fGenerateBinned, protoData, forceEvents);
       info->setCopyData(false);
   }
   return info->generate(weightVar_, protoData, forceEvents);
}
