namespace std{} 
using namespace std;
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/ExprAlgo/interface/AlgoCheck.h"

bool AlgoCheck::check(parS_type & ps, parE_type & pe, string & err)
{
  bool nok = false; // not ok - flag, negative logic inside the function!
  
  // first expression parameters
  constraintsE_type::iterator e_it = constraintsE_.begin();
  for (; e_it != constraintsE_.end(); ++e_it) {
    
    parE_type::iterator pit = pe.find(e_it->first);
    if (e_it->second.use_) { // required parameter
      
      if (pit != pe.end()) {
        nok |= checkBounds(pit,e_it,err);
      }
      else { // missing required parameter
        err += string("\tmissing required parameter: ") + e_it->first + string("\n");
	nok |= true;
      }
      	
    }
    
    else { // optional parameter
           // an optional parameter MUST be here if a default-value is specified !
      
      DCOUT('e', "AlgoCheck::check(): optional param=" << e_it->first);
      
      
      if (pit != pe.end()) { // optional parameter present 
                             // (can be a possible default-val or a user-supplied-val)
        nok |= checkBounds(pit,e_it,err);
      }	
      else { // optional parameter absent although an default value is specified. 
             // provide the default value at least minOccurs-times
        if (e_it->second.default_) {
	  for(int j=0; j<e_it->second.minOccurs_; ++j)
	    pe[e_it->first].push_back(e_it->second.defaultVal_);
	}
      }
    }
  
  }
  
  // then string parameters
  // FIXME: implement AlgoCheck::check string-valued parameter checks! 

  constraintsS_type::iterator s_it = constraintsS_.begin();
  for (; s_it != constraintsS_.end(); ++s_it) {
    
    parS_type::iterator sit = ps.find(s_it->first);
    if (s_it->second.use_) { // required parameter
      
      if (sit != ps.end()) {
	nok |= checkStrings(sit, s_it, err);
      }
      else { // missing required parameter
        err += string("\tmissing required parameter: ") + s_it->first + string("\n");
	nok |= true;
      }
    }
    else { // optional parameter absent although a default value is specified. 
           // provide the default value at least minOccurs-times

      if (s_it->second.default_) {
	for (int j = 0; j < s_it->second.minOccurs_; ++j)
	  ps[s_it->first].push_back(s_it->second.defaultVal_);
      }
    }
  }
  return !nok;
}


// checks occurences and bounds; returns true, if it is NOT OK !
bool AlgoCheck::checkBounds(parE_type::iterator pit, 
                            constraintsE_type::iterator cit, 
			    string & err 
			    )
{
   bool nok = false;
   
   DCOUT('e', "AlgoCheck::checkBounds(), par-name=" << pit->first );
   DCOUT('e', "                          minOccurs  = " << cit->second.minOccurs_ );
   DCOUT('e', "                          maxOccurs  = " << cit->second.maxOccurs_ );
   DCOUT('e', "                          par-occurs = " << pit->second.size() );
   // occurences
   if (cit->second.minOccurs_ > int(pit->second.size()) ) {
     err += "\tpar. " + cit->first 
	 +  " occurs=" + d2s(pit->second.size()) + string(" < minOccurs; minOccurs=")
	 +  d2s(cit->second.minOccurs_) + string(" maxOccurs=") 
	 +  d2s(cit->second.maxOccurs_)
	 +  "\n";
     nok |= true; 
     DCOUT('e', "                          VIOLATION of minOccurs");
   }
   //if (!nok) {
     if (cit->second.maxOccurs_ < int(pit->second.size()) ) {
     err += "\tpar. " + cit->first 
	 +  " occurs=" + d2s(pit->second.size()) + string(" > maxOccurs; minOccurs=")
	 +  d2s(cit->second.minOccurs_) + string(" maxOccurs=") 
	 +  d2s(cit->second.maxOccurs_)
	 +  "\n";
       nok |= true;   
       DCOUT('e', "                          VIOLATION of maxOccurs");
     }
   //} 
   
   // min, max
  // if (!nok) { // only check bounds if occurences were right!
     int c=0; // error count
     int i=0;
     int m=int(pit->second.size());
     for (;i<m;++i) {
       if (pit->second[i] < cit->second.min_) {
         err += "\tpar. " + cit->first + string("[") + d2s(i) + string("]=")
	     +  d2s(pit->second[i]) + string(" < min; min=")
	     +  d2s(cit->second.min_) + string(" max=") +  d2s(cit->second.max_)
	     +  "\n";
         c++; // error count
       }
       if (pit->second[i] > cit->second.max_) {
         err += "\tpar. " + cit->first + string("[") + d2s(i) + string("]=")
	     +  d2s(pit->second[i]) + string(" > max; min=")
	     +  d2s(cit->second.min_) + string(" max=") +  d2s(cit->second.max_)
	     +  "\n";
         c++; // error count       
       }
     }
     if (c)
      nok |= true;
   //}
   return nok;
}


#include <cstdio>
string AlgoCheck::d2s(double x)
{
  char buffer [25]; 
  /*int n =*/ sprintf(buffer,"%g",x);
  return string(buffer);
}

bool AlgoCheck::checkStrings(parS_type::iterator sit, 
			  constraintsS_type::iterator cit, 
			  string & err)
{
   // occurences
   bool nok = false;
   if (cit->second.minOccurs_ > int(sit->second.size()) ) {
     err += "\tpar. " + cit->first 
	 +  " occurs=" + d2s(sit->second.size()) + string(" < minOccurs; minOccurs=")
	 +  d2s(cit->second.minOccurs_) + string(" maxOccurs=") 
	 +  d2s(cit->second.maxOccurs_)
	 +  "\n";
     nok |= true; 
     DCOUT('e', "                          VIOLATION of minOccurs");
   }
   //if (!nok) {
     if (cit->second.maxOccurs_ < int(sit->second.size()) ) {
     err += "\tpar. " + cit->first 
	 +  " occurs=" + d2s(sit->second.size()) + string(" > maxOccurs; minOccurs=")
	 +  d2s(cit->second.minOccurs_) + string(" maxOccurs=") 
	 +  d2s(cit->second.maxOccurs_)
	 +  "\n";
       nok |= true;   
       DCOUT('e', "                          VIOLATION of maxOccurs");
     }
   //} 

  int i = 0;
  int m = int(sit->second.size());  
  for (; i<m;++i){

    string s = sit->second[i];
    int ss = int(s.size());

    int j = 0;
    for (; j < ss; ++j) {
      if (s[j] != ' '){
	break;
      }
      if (j == ss && s[j] == ' ') {
	err += "\tpar. " + cit->first + string("[") + d2s(i) + string("] is blank");
	nok |= true;
	DCOUT('e', "                         VIOLATION of requiring non-blank strings");
      }
    }
  }

  return nok;

}


