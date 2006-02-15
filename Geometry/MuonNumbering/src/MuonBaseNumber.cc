#include "Geometry/MuonNumbering/interface/MuonBaseNumber.h"
#include <iostream>

//#define DEBUG

void MuonBaseNumber::addBase(LevelBaseNumber num){
  basenumber_type::iterator cur=sortedBaseNumber.begin();
  basenumber_type::iterator end=sortedBaseNumber.end();

  // do a small check if level is already occupied

  while (cur!=end) {
    if (num.level()==(*cur).level()) {
      std::cout << "MuonBaseNumber::AddBase tries to add "
	   <<num.level()<<" "
	   <<num.super()<<" "
	   <<num.base()
	   <<" to existing level "
	   <<(*cur).level()<<" "
	   <<(*cur).super()<<" "
	   <<(*cur).base()
	   <<std::endl;
    }
    cur++;
  }

  cur=sortedBaseNumber.begin();
  while (cur!=end) {
    if (num.level()<(*cur).level()) break;
    cur++;
  }
  sortedBaseNumber.insert(cur,num);

#ifdef DEBUG
  cur=sortedBaseNumber.begin();
  end=sortedBaseNumber.end();
  std::cout << "MuonBaseNumber::AddBase ";
  while (cur!=end) {
    std::cout<<(*cur).level()<<" ";
    std::cout<<(*cur).super()<<" ";
    std::cout<<(*cur).base();
    std::cout<<",";
    cur++;
  }
  std::cout <<std::endl;
#endif

}

void MuonBaseNumber::addBase(const int level,const int super,const int base){
  LevelBaseNumber num(level,super,base);
  addBase(num);
}

int MuonBaseNumber::getLevels() const {
  return sortedBaseNumber.size();
}

int MuonBaseNumber::getSuperNo(int level) const {
  basenumber_type::const_iterator cur=sortedBaseNumber.begin();
  basenumber_type::const_iterator end=sortedBaseNumber.end();
  while (cur!=end) {
    if ((*cur).level()==level) {
      return (*cur).super();
    }
    cur++;
  }
  return 0;
}

int MuonBaseNumber::getBaseNo(int level) const {
  basenumber_type::const_iterator cur=sortedBaseNumber.begin();
  basenumber_type::const_iterator end=sortedBaseNumber.end();
  while (cur!=end) {
    if ((*cur).level()==level) {
      return (*cur).base();
    }
    cur++;
  }
  return 0;
}





