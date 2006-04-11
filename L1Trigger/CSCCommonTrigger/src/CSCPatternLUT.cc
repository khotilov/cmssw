#include <L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h>

int CSCPatternLUT::getBendValue(int pattern)
{
  int BendList[CSCConstants::NUM_CLCT_PATTERNS] = {0, 3, -3, 2, -2, 1, -1, 0};
  return BendList[pattern];
}

double CSCPatternLUT::getPosition(int pattern)
{
  double PositionList[CSCConstants::NUM_CLCT_PATTERNS] = {0.0, 0.0, 0.0, -0.41, 0.41, 0.42, -0.42, 0.0};
  return PositionList[pattern];
}
