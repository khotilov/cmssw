#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"


HBHEDataFrame::HBHEDataFrame() : id_(0), 
				 size_(0),
				 hcalPresamples_(0)
{
}

HBHEDataFrame::HBHEDataFrame(const HcalDetId& id) : 
  id_(id), 
  size_(0),
  hcalPresamples_(0)
{
  // TODO : test id for HcalBarrel or HcalEndcap
}

void HBHEDataFrame::setSize(int size) {
  if (size>MAXSAMPLES) size_=MAXSAMPLES;
  else if (size<=0) size_=0;
  else size_=size;
}
void HBHEDataFrame::setPresamples(int ps) {
  hcalPresamples_|=ps&0xF;
}
void HBHEDataFrame::setReadoutIds(const HcalElectronicsId& eid) {
  electronicsId_=eid;
}

bool HBHEDataFrame::validate(int firstSample, int nSamples) const {
  int capid=-1;
  bool ok=true;
  for (int i=0; ok && i<nSamples && i+firstSample<size_; i++) {
    if (data_[i+firstSample].er() || !data_[i+firstSample].dv()) ok=false;
    if (i==0) capid=data_[i+firstSample].capid();
    if (capid!=data_[i+firstSample].capid()) ok=false;
    capid=(capid+1)%4;
  }
  return ok;
}

void HBHEDataFrame::setZSInfo(bool unsuppressed, bool markAndPass) {
  if (markAndPass) hcalPresamples_|=0x10;
  if (unsuppressed) hcalPresamples_|=0x20;
}


std::ostream& operator<<(std::ostream& s, const HBHEDataFrame& digi) {
  s << digi.id() << " " << digi.size() << " samples  " << digi.presamples() << " presamples ";
  if (digi.zsUnsuppressed()) s << " zsUS ";
  if (digi.zsMarkAndPass()) s << " zsM&P ";
  s << std::endl;
  for (int i=0; i<digi.size(); i++) 
    s << "  " << digi.sample(i) << std::endl;
  return s;
}
  

