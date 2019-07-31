#ifndef MTDDETECTORDIMENSIONS_HH
#define MTDDETECTORDIMENSIONS_HH

#include <string>


namespace MTDDetectorDimensions{
  void setDimensions();
  double getDimension(std::string const& subdet);
}


#endif
