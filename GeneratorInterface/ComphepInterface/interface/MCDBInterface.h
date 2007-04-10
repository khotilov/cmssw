//************************************
//* Generators interface with the LCG MCDB
//*
//* Hector Naves Sordo 
//* 
//* First version: 25/10/06 
//* 
//************************************

#include "FWCore/Framework/interface/Event.h"
#include "Utilities/StorageFactory/interface/StorageFactory.h"
#include "Utilities/StorageFactory/interface/StorageAccount.h"
#include "FWCore/PluginManager/interface/PluginManager.h"
#include "SealBase/Storage.h"
#include "SealBase/DebugAids.h"
#include "SealBase/Signal.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace edm;
using namespace std;
using namespace seal;

// Makes a local copy of a CASTOR file.
// This code is a modified version of 
// Utilities/StorageFactory/test/any.cpp by Vincenzo Innocente
//
void mcdbGetInputFile(string  &madGraphInputFile, int &mcdbArticleID);


