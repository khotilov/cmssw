
// $Id: LumiSummary.cc,v 1.2 2008/02/06 20:57:26 wdd Exp $

#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include <iomanip>

using namespace std;

bool
LumiSummary::isProductEqual(LumiSummary const& next) {
  return (avginsdellumi_ == next.avginsdellumi_ &&
          avginsdellumierr_ == next.avginsdellumierr_ &&
          lumisecqual_ == next.lumisecqual_ &&
          deadfrac_ == next.deadfrac_ &&
          lsnumber_ == next.lsnumber_ &&
          l1ratecounter_ == next.l1ratecounter_ &&
          l1scaler_ == next.l1scaler_ &&
          hltratecounter_ == next.hltratecounter_ &&
          hltscaler_ == next.hltscaler_ &&
          hltinput_ == next.hltinput_);
}

std::ostream& operator<<(std::ostream& s, const LumiSummary& lumiSummary) {

  const std::vector<int>& l1ratecounter  = lumiSummary.l1RateCounter();
  const std::vector<int>& l1scaler       = lumiSummary.l1Scaler();
  const std::vector<int>& hltratecounter = lumiSummary.hltRateCounter();
  const std::vector<int>& hltscaler      = lumiSummary.hltScaler();
  const std::vector<int>& hltinput       = lumiSummary.hltInput();

  unsigned int maxSize = l1ratecounter.size();
  if (l1scaler.size() > maxSize) maxSize = l1scaler.size();
  if (hltratecounter.size() > maxSize) maxSize = hltratecounter.size();
  if (hltscaler.size() > maxSize) maxSize = hltscaler.size();
  if (hltinput.size() > maxSize) maxSize = hltinput.size();

  s << "\nDumping LumiSummary\n\n";

  s << "  avgInsDelLumi = " << lumiSummary.avgInsDelLumi() << "\n";
  s << "  avgInsDelLumiErr = " << lumiSummary.avgInsDelLumiErr() << "\n";
  s << "  lumiSecQual = " << lumiSummary.lumiSecQual() << "\n";
  s << "  deadFrac = " << lumiSummary.deadFrac() << "\n";
  s << "  liveFrac = " << lumiSummary.liveFrac() << "\n";
  s << "  lsNumber = " << lumiSummary.lsNumber() << "\n";
  s << "  avgInsRecLumi = " << lumiSummary.avgInsRecLumi() << "\n";
  s << "  avgInsRecLumiErr = " << lumiSummary.avgInsRecLumiErr() << "\n\n";

  s << setw(15) << "l1ratecounter";
  s << setw(15) << "l1scaler";
  s << setw(15) << "hltratecounter";
  s << setw(15) << "hltscaler";
  s << setw(15) << "hltinput";
  s << "\n";

  for (unsigned int i = 0; i < maxSize; ++i) {

    s << setw(15);
    i < l1ratecounter.size() ? s << l1ratecounter[i] : s << " ";

    s << setw(15);
    i < l1scaler.size() ? s << l1scaler[i] : s << " ";

    s << setw(15);
    i < hltratecounter.size() ? s << hltratecounter[i] : s << " ";

    s << setw(15);
    i < hltscaler.size() ? s << hltscaler[i] : s << " ";

    s << setw(15);
    i < hltinput.size() ? s << hltinput[i] : s << " ";

    s << "\n";
  }
  return s << "\n";
}
