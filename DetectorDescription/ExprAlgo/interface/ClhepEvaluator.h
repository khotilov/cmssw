#ifndef ClhepEvaluator_h
#define ClhepEvaluator_h


#include <vector>
#include <string>
#include "DetectorDescription/ExprAlgo/interface/ExprEvalInterface.h"
#include "CLHEP/Evaluator/Evaluator.h"


class ClhepEvaluator : public ExprEvalInterface
{
public:
 ClhepEvaluator();
  
  ~ClhepEvaluator();
  
  void set(const string & ns, const string & name, const string & exprValue);
   
  double eval(const string & ns, const string & expr);
  
  bool isDefined(const string & ns, //< current namespace
                 const string & name //< name of the variable inside current namespace
		 );

  //! access to the clhep-implementation of the dictionary variables
  const std::vector<std::string> & variables() const { return variables_;}
  const std::vector<std::string> & values() const { return values_;}
  
  //! evaluations using directly the CLHEP-evaluator
  /** expression must be an expression compatible with the CLHEP-Evaluator syntax */
  double eval(const char * expression);
  
  //! filling the clhep-implementation of the dictionary
  void set(const string & name, const string & value);
  
  void clear();

private:
  void prepare(const std::string & ns,         // input
               const std::string & name,       // input
	       const std::string & exprValue,  // input
               std::string & nameResult,       // output
	       std::string & valResult) const; // output
	       
  void throwex(const std::string & ns, 
               const std::string & name, 
	       const std::string & expr,
	       const std::string & reason,
	       int idx=0) const;
	       
  void checkname(const std::string & name) const;  	       
  
  HepTool::Evaluator evaluator_;  	
  std::vector<std::string> variables_;
  std::vector<std::string> values_;         
};

#endif
