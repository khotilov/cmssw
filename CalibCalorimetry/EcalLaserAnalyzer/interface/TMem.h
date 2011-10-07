#ifndef TMem_H
#define TMem_H

#include<vector>
using namespace std;

class TMem 
{

 private:

  int _fedid;
  vector <int> _memFromDcc;
  
  void init(int);

 public:


  // Default Constructor, mainly for Root
  TMem();

  // Constructor
  TMem(int);

  // Destructor: Does nothing
  virtual ~TMem();

  bool isMemRelevant(int);
  int Mem(int, int);

  //ClassDef(TMem,1)
};

#endif
