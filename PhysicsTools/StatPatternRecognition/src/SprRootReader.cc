// $Id: SprRootReader.cc,v 1.9 2007/10/10 20:21:31 narsky Exp $

#include "PhysicsTools/StatPatternRecognition/interface/SprExperiment.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprRootReader.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprData.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprAbsFilter.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprEmptyFilter.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprPreFilter.hh"

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TObjArray.h>

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <memory>

using namespace std;

SprRootReader::SprRootReader(SprPreFilter* filter)
  : 
  SprAbsReader(filter),
  treeName_(),
  leafNames_(),
  fileObjects_(),
  hasSpecialClassifier_(false),
  classifierVarName_()
{}

// parses the text file to read names of root files
// defers reading of those to readRootObjects()
SprAbsFilter* SprRootReader::read(const char* filename)
{
    ifstream file(filename);
    if (not file) {
        cerr << "Unable to open " << filename << endl;
        return 0;
    }

    string line;
    double weight = 1.0;
    // if the weight is never set, we can save some time
    bool weightHasChanged = false;
    cout << "Parsing File: " << filename << endl;
    while (getline(file, line)) {
        if (line.find('#') != string::npos) {
            line.erase(line.find_first_of('#'));
        }
        if (line.find_first_not_of(' ') == string::npos)
            continue;

        istringstream inString(line);
        vector<string> lineFields;
        string fieldDummy;
        while (inString >> fieldDummy)
            lineFields.push_back(fieldDummy);

	assert( lineFields.size() > 1 );

        if (lineFields.at(0) == "Tree:") {
            assert(treeName_.length() == 0);
            treeName_ = lineFields.at(1);
	}else if (lineFields.at(0) == "ClassVariable:") {
            //Accept variable name of TrueClass	
	    if(hasSpecialClassifier_){
  	      cout<<"WARNING - True class variable was already chosen as "
		  <<classifierVarName_<<" will be overwritten to "
		  <<lineFields.at(1)
		  <<"\nPlease change your Run File"<<endl;
	    }	  	  
            hasSpecialClassifier_ = true;
            classifierVarName_ = lineFields.at(1);
        } else if (lineFields.at(0) == "WeightVariable:") {
            //Accept variable name of TrueClass
            weightHasChanged = true;
            for (int dummy=1; dummy<lineFields.size(); ++dummy) {
                weightLeafNames_.push_back(lineFields.at(dummy));
            }

            if(weightLeafNames_.size() == 0){
              cout<<"WARNING - Parsing error concerning weight variables"
                  <<"\nPlease change your Run File"<<endl;
            }
        } else if (lineFields.at(0) == "Leaves:") {
            for (int dummy=1; dummy<lineFields.size(); ++dummy) {
                leafNames_.push_back(lineFields.at(dummy));
            }
        } else if (lineFields.at(0) == "Weight:") {
            weightHasChanged = true;
            istringstream s(lineFields.at(1));
            s >> weight;
        } else if (lineFields.at(0) == "File:") {

  	    assert( lineFields.size() > 1 );

            FileInfo thisFile;
            thisFile.name = lineFields.at(1);

	    thisFile.start = 0;
	    thisFile.end = -1;

	    if( lineFields.size() > 2 ) {
	      istringstream 
		dummyIn(string(lineFields.at(2),
			       0,
			       lineFields.at(2).find_first_of('-')));
	      if( !(dummyIn >> thisFile.start) ) {
                thisFile.start = 0;
	      }
	      dummyIn.clear();
	      dummyIn.str(string(lineFields.at(2), 
				 lineFields.at(2).find_first_of('-')+1, 
				 string::npos));
	      if( !(dummyIn >> thisFile.end) ) {
                thisFile.end = -1;
	      }
	      dummyIn.clear();
	    }

	    thisFile.fileClass = 0;

	    if( lineFields.size() > 3 ) {
	      istringstream dummyIn(lineFields.at(3));
	      if (not (dummyIn >> thisFile.fileClass)) {
                thisFile.fileClass = 0;
                cout << dummyIn.get();
	      }
	    }

	    thisFile.weight = weight;
            fileObjects_.push_back(thisFile);

            cout << "For tree " << treeName_
		 << " found file: " << thisFile.name
		 << " start: " << thisFile.start
		 << " end: " << thisFile.end
		 << " class: " << thisFile.fileClass
		 << " weight: " << thisFile.weight
		 << endl;
        }
    }

    if(hasSpecialClassifier_){
      cout << "True class value is given by leaf " 
	   << classifierVarName_ << endl;
    }

    if(weightLeafNames_.size()){
        cout<<"A variable determined weight has been chosen, the value"
            <<" assigned to ";
	for(int i = 0; i < weightLeafNames_.size(); i++){
           if(i%5 == 0) cout<<"\n\t";
	   if(i == 0)  cout<<weightLeafNames_[i];
           else cout<<" * "<<weightLeafNames_[i];
        }
	cout<<"\n will be used for the weight."<<endl;
    }
    
    return readRootObjects(weightHasChanged);
}

SprAbsFilter* SprRootReader::readRootObjects(bool needToCalcWeights)
{
    auto_ptr<SprData> data(new SprData);
    vector<double> weights;

    // set up pre-filter vars
    if( filter_!=0 && !filter_->setVars(leafNames_) ) {
      cerr << "Unable to apply pre-filter requirements." << endl;
      return 0;
    }

    // get a new list of variables
    vector<string> transformed;
    if( filter_ != 0 ) {
      if( !filter_->transformVars(leafNames_,transformed) ) {
	cerr << "Pre-filter is unable to transform variables." << endl;
	return 0;
      }
    }
    if( transformed.empty() ) transformed = leafNames_; 

    // set up data vars
    if( !data->setVars(transformed) ) {
      cerr << "Unable to set variable list for input data." << endl;
      return 0;
    }

    // fill out data
    for(vector<FileInfo>::iterator fileIter = fileObjects_.begin()
            ; fileIter != fileObjects_.end()
            ; ++fileIter) {
        cout << "Reading File: " << fileIter->name;
        TFile f(fileIter->name.c_str());
        TTree* tree = dynamic_cast<TTree*>(f.Get(treeName_.c_str()));
        if(tree == 0) {
           cout<<"Tree "<<treeName_<<" not found in file "
	       <<fileIter->name.c_str()<<endl;
           continue;
        }
	if (fileIter->end < 0)
            fileIter->end = tree->GetEntries();
        cout << " (" << fileIter->end - fileIter->start << " events)" << endl;
        map<string, TLeaf*> leaves;
        for (vector<string>::const_iterator leafIter = leafNames_.begin()
                ; leafIter != leafNames_.end()
                ; ++leafIter) {

            TLeaf* tempLeaf = tree->GetLeaf(leafIter->c_str());

            if(tempLeaf == 0){
                 cerr<<"No Leaf associated with variable "
                     <<leafIter->c_str()<<" - probably a typo"<<endl;
                 abort();
            }
            leaves.insert(make_pair(*leafIter, tempLeaf));

        }
        for (int iEvent=fileIter->start; iEvent<fileIter->end; ++iEvent) {
            tree->GetEntry(iEvent);
            vector<double> row;
            for (vector<string>::const_iterator leafIter = leafNames_.begin()
                 ; leafIter != leafNames_.end()
                 ; ++leafIter) {
                    // Take always the first entry
                    row.push_back(leaves[*leafIter]->GetValue(0));
            }

            int assignedClass = fileIter->fileClass;	    
            if(hasSpecialClassifier_){
		TLeaf* classLeaf = tree->GetLeaf(classifierVarName_.c_str());
		if(classLeaf == 0){
		  cerr<<"No Leaf associated with classifier - probably"
		      <<" a typo - will assign classification 0"<<endl;
		  abort();
		}else{
		  assignedClass = (int) classLeaf->GetValue(0);
                }
	    }
    
            float assignedWeight = fileIter->weight;
            for (int i = 0; i < weightLeafNames_.size(); i++)
            {   
                TLeaf* tempLeaf = tree->GetLeaf(weightLeafNames_[i].c_str());
		if(tempLeaf == 0){
                  cerr<<"No Leaf associated with variable "
		      <<weightLeafNames_[i]<<" - probably"
                      <<" a typo - please fix this"<<endl;
                  abort();
                }else{
                  assignedWeight *= (float) tempLeaf->GetValue(0);
                } 
            }

	    if( filter_!=0 && !filter_->pass(assignedClass,row) ) continue;
	    if( filter_ != 0 ) {
	      vector<double> vNew;
	      if( filter_->transformCoords(row,vNew) ) {
		data->insert(assignedClass,vNew);
		weights.push_back(assignedWeight);
		continue;
	      }
	      cerr << "Pre-filter is unable to transform coordinates." << endl;
	      return 0;
	    }
	    data->insert(assignedClass,row);
            weights.push_back(assignedWeight);
        }
        f.Close();
    }
    if (needToCalcWeights)
        return new SprEmptyFilter(data.release(), weights, true);

    return new SprEmptyFilter(data.release(), true);
}


bool SprRootReader::chooseVars(const std::set<std::string>& vars)
{
  cerr << "Unable to choose variables: "
       << "SprRootReader::chooseVars() not implemented." << endl;
  return false;
}


bool SprRootReader::chooseAllBut(const std::set<std::string>& vars)
{
  cerr << "Unable to choose variables: "
       << "SprRootReader::chooseAllBut() not implemented." << endl;
  return false;
}


void SprRootReader::chooseAll()
{
  vector<FileInfo>::iterator fileIter = fileObjects_.begin();
  TFile f(fileIter->name.c_str());
  TTree* tree = dynamic_cast<TTree*>(f.Get(treeName_.c_str()));
  
  if(tree == 0) {
    cerr<<"Tree "<<treeName_<<" not found in file "
	<<fileIter->name.c_str()<<endl;
    cerr<<"No variables will be selected."<<endl;
    return;
  }
  
  TObjArray* leafArray = tree->GetListOfLeaves();
  TIter leafIter(leafArray);
  leafIter.Reset();
  
  leafNames_.clear();
  TLeaf* thisLeaf = 0;
  while( (thisLeaf = (TLeaf*)leafIter.Next()) != 0 ){
    leafNames_.push_back(thisLeaf->GetName());
  }
}
