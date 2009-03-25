
#include <math.h>

#include "PtrelSolver.h"

#include "TClass.h"
#include "TDirectory.h"
#include "TFractionFitter.h"
#include "TKey.h"
#include "TMinuit.h"
#include "TObjString.h"

#include "PtrelUtils.h"


ClassImp(PtrelSolver)


char const * PtrelSolver::directory = "/Histograms/muon_in_jet";


void PtrelSolver::templates(const char * filename)
{
    if (templates_) templates_->Close();

    GetSafely(templates_, TFile::Open(filename, "READ"))

    char name[256];

    for (Int_t i = 0; i < Dependency::Dimension; ++i)
    {
        sprintf(name, "/parameters/%s", Dependency::Name[i]);
        TObjString * value = (TObjString *) templates_->Get(name);
        if (value)
        {
            rebin_[i] = value->String().Atoi();
            Info(__FUNCTION__, "Loading rebinning factor %d for %s", rebin_[i], Dependency::Name[i]);
        }
    }
}


bool PtrelSolver::measure(
    TFile * input,
    TFile * output,
    const char * name,
    ValueVector & vVector,
    ValueVector & eVector
)
{
    CreateSafelyZero(TH2D, histogram, readTH2(input, name))
    CallSafelyZero(measure(output, histogram, vVector, eVector))
    return true;
}


bool PtrelSolver::measure(
    TFile * input,
    TFile * output,
    TPRegexp pattern,
    StringVector & hVector,
    ValueMatrix & vMatrix,
    ValueMatrix & eMatrix
)
{
    // Return status
    bool status = false;

    // Clean the containers
    hVector.clear();
    vMatrix.clear();
    eMatrix.clear();

    // Move to the directory with 2d histograms
    input->cd(directory);

    // Loop over all keys in this directory
    TIter nextkey( gDirectory->GetListOfKeys() );

    TKey * key;

    //char name[256];

    while (( key = (TKey*)nextkey() ))
    {
        // Measure those 2D histogram with name fit the pattern
        TObject * object = key->ReadObj();
        if ( object->IsA()->InheritsFrom( "TH2" ) )
            if ( TString(object->GetName() ).Contains(pattern) )
            {
                // Update status
                status = true;

                // Temporal container for efficiencies
                TH2 * histogram;
                ValueVector values, errors;

                // Process 2D histogram associated to the keyword
                GetSafelyZero(histogram, processTH2(object))

                // Measuring efficiencies
                if (!measure(output, histogram, values, errors))
                {
                    Warning(__FUNCTION__, "Measure failed skipping it");
                    continue;
                }

                // Appending the results
                hVector.push_back(TString(histogram->GetName()));
                vMatrix.push_back(values);
                eMatrix.push_back(errors);
            }
    };

    if (!status) Error(__FUNCTION__, "Non matching histograms were found");

    return status;
}


bool PtrelSolver::measure(
    TFile * output,
    TH2 * histogram2D,
    ValueVector & vVector,
    ValueVector & eVector
)
{
    char name[256];

    // Creating sub directory for templates
    if (!output->FindKey("fits")) output->mkdir("fits");

    // Clearing the containers
    vVector.clear();
    eVector.clear();

    // Number of bins
    Int_t nbins = histogram2D->GetNbinsX();

    // Loop over the bins
    for (Int_t i = 1; i <= nbins; ++i)
    {
        // Information
        Info(__FUNCTION__, "Measuring flavor content using fit option %s in bin %d", Fit::Label[fittype_], i);

        if (fittype_ == Fit::histograms)
        {
            // Get the proper set of histograms
            CreateSafelyZero(TObjArray, templates, combinedHistograms(histogram2D->GetName(), i))

            // Set histogram name by the formula = data_x_bin (transient)
            sprintf(name, "data_%s_%d", histogram2D->GetName(), i);
            // Project histogram
            TH1D * histogram1D = histogram2D->ProjectionY(name, i, i, "e");

            // Temporal container with the results
            TVectorD values, errors;

            // Fitting the histogram
            CallSafelyZero(fit(output, histogram1D, templates, values, errors))

            // Collecting the results
            vVector.push_back(values);
            eVector.push_back(errors);
        }
        else
        {
            // Get the proper combined function
            CreateSafelyZero(TF1, templates, combinedFunction(histogram2D->GetName(), i))

            // Set histogram name by the formula = fit_x_bin
            sprintf(name, "fit_%s_%d", histogram2D->GetName(), i);
            // Project histogram
            TH1D * histogram1D = histogram2D->ProjectionY(name, i, i, "e");
            // Basic setup
            ptrelHistogramSetup(histogram1D);

            // Temporal container with the results
            TVectorD values, errors;

            // Fitting the histogram
            CallSafelyZero(fit(histogram1D, templates, values, errors))

            // Collecting the results
            vVector.push_back(values);
            eVector.push_back(errors);

            // Saving the fitting.
            output->cd();
            output->cd("fits");
            histogram1D->Write();
        }
    }
    return true;
}


TH2 * PtrelSolver::readTH2(TFile * file, const char * name) const
{
    CreateSafelyZero(TObject, object, file->Get(name))
    CreateSafelyZero(TH2, histogram, processTH2(object))
    return histogram;
}


TH2 * PtrelSolver::processTH2(TObject * object) const
{
    // Cast the object pointer into 2D histogram
    TH2 * histogram = (TH2*) object;

    // Rebinning in ptrel
    if ( rebin_[Dependency::ptrel] != 1 )
    {
        Info(__FUNCTION__, "Rebinning ptrel a factor %d in %s", rebin_[Dependency::ptrel], histogram->GetName());
        histogram->RebinY(rebin_[Dependency::ptrel]);
    }

    // Check the histogram dependency and rebin
    for (Int_t i = 1; i < Dependency::Dimension; ++i)
        if ( containsIdentifier(histogram->GetName(),Dependency::Name[i]) && rebin_[i] != 1)
        {
            Info(__FUNCTION__, "Rebinning %s a factor %d in %s", Dependency::Name[i], rebin_[i], histogram->GetName());
            histogram->RebinX(rebin_[i]);
        }

    return histogram;
}


bool PtrelSolver::fit(TH1 * histogram, TF1 * function, TVectorD & values, TVectorD & errors)
{
    // Check the number of flavor is not lower than two, otherwise cancel the fitting.
    if ( fitFlavors_.size() < 2 )
    {
        Error(__FUNCTION__, "Number of flavor to fit is lower than two.");
        return false;
    }

    // Check for initial values (otherwise using default)
    if (values.GetNoElements() != (Int_t) fitFlavors_.size())
    {
        Warning(__FUNCTION__, "Inconsistent dimension for initial values using default");
        values.ResizeTo(fitFlavors_.size());
        values.Zero();
        for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
            values(i) = 1./fitFlavors_.size();
    }

    // Numerical integration of data histogram in the interval of the function
    Double_t integral = histogram->Integral(
                            histogram->FindBin(function->GetXmin()),
                            histogram->FindBin(function->GetXmax()),
                            "width"
                        );

    // Calculation of normalization parameters
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
    {
        function->SetParameter(i, integral * values(i));
        Info(__FUNCTION__, "Fitting with initial fraction %f for %s", values(i), Flavor::Name[fitFlavors_[i]]);
    }

    // Fitting combined function against data
    histogram->Fit(function, "Q", "", function->GetXmin(), function->GetXmax());

    // Information
    Info(__FUNCTION__, "Fitting %s chi2/ndf = (%f/%d)", histogram->GetName(), function->GetChisquare(), function->GetNDF());

    // Temporal containers for the results
    TVectorD values_(fitFlavors_.size());
    TMatrixD covariance_(fitFlavors_.size(), fitFlavors_.size());

    // Getting the results for the fitting
    function->GetParameters(values_.GetMatrixArray());
    // Get the covariance matrix of parameter errors
    gMinuit->mnemat(covariance_.GetMatrixArray(), fitFlavors_.size());

    // Copy semantics for the values
    values = values_;

    // Total number of event including correction because of weights
    Double_t numberEvents = histogram->Integral();

    // Sum of the x values
    Double_t sumx = values.Sum();

    // Derivative matrix for error computation
    TMatrixD derivative(fitFlavors_.size(), fitFlavors_.size());
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
        for (Int_t j = 0; j < (Int_t) fitFlavors_.size(); ++j)
        {
            if (i == j)
                derivative(j,i) = numberEvents/sumx - numberEvents*values(i)/(sumx*sumx);
            else
                derivative(j,i) = - numberEvents*values(i)/(sumx*sumx);
        }

    // Scaling to the number of events
    values *= numberEvents/sumx;

    // Error computation
    TMatrixD temporal (covariance_, TMatrixD::kMult, derivative);
    TMatrixD covariance (derivative.T(), TMatrixD::kMult, temporal);

    errors.ResizeTo(fitFlavors_.size());
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
        errors = sqrt(covariance(i,i));

    return true;
}


bool PtrelSolver::fit(TFile * output, TH1 * histogram, TObjArray * templates, TVectorD & values, TVectorD & errors)
{
    // Check the number of flavor is not lower than two, otherwise cancel the fitting.
    if ( fitFlavors_.size() < 2 )
    {
        Error(__FUNCTION__, "Number of flavor to fit is lower than two.");
        return false;
    }

    // Check if data histogram is empty
    if (!histogram->GetEntries())
    {
        Warning(__FUNCTION__, "Empty data histogram !");
        return false;
    }

    // Creating a fraction fitter object
    TFractionFitter fit(histogram, templates);

    // Set the constrains for all fraction values.
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
        fit.Constrain(i, -1.0, 2.0);

    // Running the fitter
    Int_t status = fit.Fit();

    // Resizing error vector
    values.ResizeTo(fitFlavors_.size());
    errors.ResizeTo(fitFlavors_.size());

    // Warning if the fit do not converge
    if (status)
    {
        Warning(__FUNCTION__, "Fitting do not converge (status code %d)", status);
        values.Zero();
        errors.Zero();
        return true;
    }

    // Information
    Info(__FUNCTION__, "Fitting %s chi2/ndf = (%f/%d)", histogram->GetName(), fit.GetChisquare(), fit.GetNDF());

    // Resizing error vector
    values.ResizeTo(fitFlavors_.size());
    errors.ResizeTo(fitFlavors_.size());

    // Total number of event including correction because of weights
    Double_t numberEvents = histogram->Integral();

    // Getting the results and scaling up to the number of events
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
    {
        Double_t value, error;
        fit.GetResult(i, value, error);
        values(i) = numberEvents * value;
        errors(i) = numberEvents * error;
    }

    // Saving the fitting.
    output->cd();
    output->cd("fits");
    TH1 * result = fit.GetPlot();
    ptrelHistogramSetup(result);
    result->Write();

    return true;
}


TF1 * PtrelSolver::combinedFunction(char const * keyword, Int_t bin)
{
    char name[256];
    TString form;
    Double_t xmin = 0;
    Double_t xmax = 0;

    // Look over the flavor functions and creating
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
    {
        // HACK to solve the problem of naming conventio
        TPRegexp p("[A-Z]{3,}$");
        if (TString(keyword).Contains(p))
        {
            TPRegexp t("[A-Z]{3,}");
            Int_t inx = TString(keyword).Index(t);
            std::string tmp(keyword);
            std::string tag(tmp.substr(inx, tmp.size()));
            std::string core(tmp.substr(0, inx-1));
            // Read the function from the file
            sprintf(name, "/functions/function_%s_%s_%s_%d", core.c_str(), Flavor::Name[fitFlavors_[i]], tag.c_str(), bin);
        }  // END HACK
        else
            sprintf(name, "/functions/function_%s_%s_%d", keyword, Flavor::Name[fitFlavors_[i]], bin);

        Info(__FUNCTION__, "Loading %s", name);
        CreateSafelyZero(TF1, function, templates_->Get(name))

        // Collect minimal and maximal values
        xmin = function->GetXmin();
        xmax = function->GetXmax();

        // Calculate the integral
        Double_t scale = 1./function->Integral(xmin, xmax);

        // Create the funcion form for combined function
        if (i != 0) form += '+';
        form += '[';
        form += i;
        form += "]*";
        form += scale;
        form += "*(";
        form += function->GetExpFormula("p");
        form += ')';
    }

    // Information
    sprintf(name, "combined_%s_%d", keyword, bin);
    Info(__FUNCTION__, "Creating combined function %s defined as %s", name, form.Data());

    // Creating function
    return new TF1(name, form.Data(), xmin, xmax);
}


TObjArray * PtrelSolver::combinedHistograms(char const * keyword, Int_t bin)
{
    char name[256];
    TString form;

    TObjArray * histograms = new TObjArray(fitFlavors_.size());

    // Look over the flavor functions and creating
    for (Int_t i = 0; i < (Int_t) fitFlavors_.size(); ++i)
    {
        // HACK to solve the problem of naming convention
        TPRegexp p("[A-Z]{3,}$");
        if (TString(keyword).Contains(p))
        {
            TPRegexp t("[A-Z]{3,}");
            Int_t inx = TString(keyword).Index(t);
            std::string tmp(keyword);
            std::string tag(tmp.substr(inx, tmp.size()));
            std::string core(tmp.substr(0, inx-1));
            // Read the function from the file
            sprintf(name, "/templates/template_%s_%s_%s_%d", core.c_str(), Flavor::Name[fitFlavors_[i]], tag.c_str(), bin);
        }  // END HACK
        else
            sprintf(name, "/templates/template_%s_%s_%d", keyword, Flavor::Name[fitFlavors_[i]], bin);

        Info(__FUNCTION__, "Loading %s", name);
        CreateSafelyZero(TH1, histogram, templates_->Get(name))
        histograms->Add(histogram);
    }

    // Return the collection of histograms
    return histograms;
}

// Use the 'flavor' template when doing Ptrel fit
void PtrelSolver::setFitFlavor(Flavor::Type flavor)
{
    fitFlavors_.push_back(flavor);
}
