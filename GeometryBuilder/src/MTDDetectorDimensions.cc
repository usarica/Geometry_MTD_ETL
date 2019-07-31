#include <cassert>
#include <map>
#include <algorithm>
#include <utility>
#include "MTDDetectorDimensions.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"


// Could use unordered_map as well
typedef std::map<std::string, double> dim_map_t;
typedef std::map<std::string, double>::iterator dim_map_it_t;
typedef std::map<std::string, double>::const_iterator dim_map_cit_t;


namespace MTDDetectorDimensions{
  dim_map_t dimension_map;

  void setETLOneSensorModuleDimensions();
  void setETLTwoSensorModuleDimensions();
  void setETLSensorServiceHybridDimensions(int const& nSensorsPerSide); // nSensorsPerSide = 3 or 6

  void printDimensions();
}


using namespace std;
using namespace CLHEP;


double MTDDetectorDimensions::getDimension(std::string const& subdet){
  dim_map_cit_t it = dimension_map.find(subdet);
  if (it != dimension_map.cend()) return it->second;
  else{
    G4cout << "MTDDetectorDimensions::getDimension: Could not find the dimension " << subdet << G4endl;
    assert(0);
    return -1;
  }
}
void MTDDetectorDimensions::setDimensions(){
  setETLOneSensorModuleDimensions();
  setETLTwoSensorModuleDimensions();
  setETLSensorServiceHybridDimensions(3);
  setETLSensorServiceHybridDimensions(6);

  printDimensions();
}

void MTDDetectorDimensions::printDimensions(){ for (auto it:dimension_map) G4cout << "Dimension " << it.first << " = " << it.second << G4endl; }

void MTDDetectorDimensions::setETLOneSensorModuleDimensions(){
  string const detbase = "ETLOneSensorModule";
  string detname;

  // AlN base plate
  detname = detbase + "_BasePlate";
  double baseplateSize_X = 28.25*mm;
  double baseplateSize_Y = 43.1*mm;
  double baseplateSize_Z = 0.79*mm;
  dimension_map[detname+"_X"] = baseplateSize_X;
  dimension_map[detname+"_Y"] = baseplateSize_Y;
  dimension_map[detname+"_Z"] = baseplateSize_Z;


  // Thermal pad
  detname = detbase + "_BaseFilm";
  double basefilmSize_X = baseplateSize_X;
  double basefilmSize_Y = baseplateSize_Y;
  double basefilmSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = basefilmSize_X;
  dimension_map[detname+"_Y"] = basefilmSize_Y;
  dimension_map[detname+"_Z"] = basefilmSize_Z;

  // ETROC
  detname = detbase + "_ETROC";
  double etrocSize_X = 22.3*mm;
  double etrocSize_Y = 20.8*mm;
  double etrocSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = etrocSize_X;
  dimension_map[detname+"_Y"] = etrocSize_Y;
  dimension_map[detname+"_Z"] = etrocSize_Z;

  // Laird film
  detname = detbase + "_LairdFilm";
  double lairdfilmSize_X = etrocSize_X;
  double lairdfilmSize_Y = etrocSize_Y*2.;
  double lairdfilmSize_Z = 0.08*mm;
  dimension_map[detname+"_X"] = lairdfilmSize_X;
  dimension_map[detname+"_Y"] = lairdfilmSize_Y;
  dimension_map[detname+"_Z"] = lairdfilmSize_Z;

  // LGAD sensor
  detname = detbase + "_LGAD";
  double lgadSize_X = 21.2*mm;
  double lgadSize_Y = 42.0*mm;
  double lgadSize_Z = 0.3*mm;
  dimension_map[detname+"_X"] = lgadSize_X;
  dimension_map[detname+"_Y"] = lgadSize_Y;
  dimension_map[detname+"_Z"] = lgadSize_Z;

  // Solder bumps
  detname = detbase + "_SolderBumps";
  double bumpsSize_X = lgadSize_X;
  double bumpsSize_Y = etrocSize_Y*2.;
  double bumpsSize_Z = 0.03*mm;
  dimension_map[detname+"_X"] = bumpsSize_X;
  dimension_map[detname+"_Y"] = bumpsSize_Y;
  dimension_map[detname+"_Z"] = bumpsSize_Z;

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  double epoxySize_X = lgadSize_X;
  double epoxySize_Y = lgadSize_Y;
  double epoxySize_Z = 0.08*mm;
  dimension_map[detname+"_X"] = epoxySize_X;
  dimension_map[detname+"_Y"] = epoxySize_Y;
  dimension_map[detname+"_Z"] = epoxySize_Z;

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  double coverplateSize_X = lgadSize_X; // Slightly wrong, should be slightly longer
  double coverplateSize_Y = lgadSize_Y;
  double coverplateSize_Z = 0.51*mm;
  dimension_map[detname+"_X"] = coverplateSize_X;
  dimension_map[detname+"_Y"] = coverplateSize_Y;
  dimension_map[detname+"_Z"] = coverplateSize_Z;

  // All of the module
  detname = detbase;
  double moduleSize_X = baseplateSize_X;
  double moduleSize_Y = baseplateSize_Y;
  double moduleSize_Z = basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+epoxySize_Z+lgadSize_Z+coverplateSize_Z;
  dimension_map[detname+"_X"] = moduleSize_X;
  dimension_map[detname+"_Y"] = moduleSize_Y;
  dimension_map[detname+"_Z"] = moduleSize_Z;
}
void MTDDetectorDimensions::setETLTwoSensorModuleDimensions(){
  string const detbase = "ETLTwoSensorModule";
  string detname;

  // AlN base plate
  detname = detbase + "_BasePlate";
  double baseplateSize_X = 56.5*mm;
  double baseplateSize_Y = 43.1*mm;
  double baseplateSize_Z = 0.79*mm;
  dimension_map[detname+"_X"] = baseplateSize_X;
  dimension_map[detname+"_Y"] = baseplateSize_Y;
  dimension_map[detname+"_Z"] = baseplateSize_Z;


  // Thermal pad
  detname = detbase + "_BaseFilm";
  double basefilmSize_X = baseplateSize_X;
  double basefilmSize_Y = baseplateSize_Y;
  double basefilmSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = basefilmSize_X;
  dimension_map[detname+"_Y"] = basefilmSize_Y;
  dimension_map[detname+"_Z"] = basefilmSize_Z;

  // ETROC
  detname = detbase + "_ETROC";
  double etrocSize_X = 22.3*mm;
  double etrocSize_Y = 20.8*mm;
  double etrocSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = etrocSize_X;
  dimension_map[detname+"_Y"] = etrocSize_Y;
  dimension_map[detname+"_Z"] = etrocSize_Z;

  // Laird film
  detname = detbase + "_LairdFilm";
  double lairdfilmSize_X = etrocSize_X*2.;
  double lairdfilmSize_Y = etrocSize_Y*2.;
  double lairdfilmSize_Z = 0.08*mm;
  dimension_map[detname+"_X"] = lairdfilmSize_X;
  dimension_map[detname+"_Y"] = lairdfilmSize_Y;
  dimension_map[detname+"_Z"] = lairdfilmSize_Z;

  // LGAD sensor
  detname = detbase + "_LGAD";
  double lgadSize_X = 21.2*mm;
  double lgadSize_Y = 42.0*mm;
  double lgadSize_Z = 0.3*mm;
  dimension_map[detname+"_X"] = lgadSize_X;
  dimension_map[detname+"_Y"] = lgadSize_Y;
  dimension_map[detname+"_Z"] = lgadSize_Z;

  // Solder bumps
  detname = detbase + "_SolderBumps";
  double bumpsSize_X = lgadSize_X;
  double bumpsSize_Y = etrocSize_Y*2.;
  double bumpsSize_Z = 0.03*mm;
  dimension_map[detname+"_X"] = bumpsSize_X;
  dimension_map[detname+"_Y"] = bumpsSize_Y;
  dimension_map[detname+"_Z"] = bumpsSize_Z;

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  double epoxySize_X = lgadSize_X;
  double epoxySize_Y = lgadSize_Y;
  double epoxySize_Z = 0.08*mm;
  dimension_map[detname+"_X"] = epoxySize_X;
  dimension_map[detname+"_Y"] = epoxySize_Y;
  dimension_map[detname+"_Z"] = epoxySize_Z;

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  double coverplateSize_X = lgadSize_X; // Slightly wrong, should be slightly longer
  double coverplateSize_Y = lgadSize_Y;
  double coverplateSize_Z = 0.51*mm;
  dimension_map[detname+"_X"] = coverplateSize_X;
  dimension_map[detname+"_Y"] = coverplateSize_Y;
  dimension_map[detname+"_Z"] = coverplateSize_Z;

  // All of the module
  detname = detbase;
  double moduleSize_X = baseplateSize_X;
  double moduleSize_Y = baseplateSize_Y;
  double moduleSize_Z = basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+epoxySize_Z+lgadSize_Z+coverplateSize_Z;
  dimension_map[detname+"_X"] = moduleSize_X;
  dimension_map[detname+"_Y"] = moduleSize_Y;
  dimension_map[detname+"_Z"] = moduleSize_Z;
}
void MTDDetectorDimensions::setETLSensorServiceHybridDimensions(int const& nSensorsPerSide){
  string const detbase = "ETL" + std::to_string(2*nSensorsPerSide) + "SensorServiceHybrid";
  string detname;

  double servicehybridSize_X = 34*mm;
  double servicehybridSize_Y = 43.1*((double) nSensorsPerSide)*mm;

  // Thermal pad
  detname = detbase + "_ThermalPad";
  double thermalpadSize_X = servicehybridSize_X;
  double thermalpadSize_Y = servicehybridSize_Y;
  double thermalpadSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = thermalpadSize_X;
  dimension_map[detname+"_Y"] = thermalpadSize_Y;
  dimension_map[detname+"_Z"] = thermalpadSize_Z;

  // Readout board
  detname = detbase + "_ReadoutBoard";
  double readoutboardSize_X = servicehybridSize_X;
  double readoutboardSize_Y = servicehybridSize_Y;
  double readoutboardSize_Z = 1.*mm;
  dimension_map[detname+"_X"] = readoutboardSize_X;
  dimension_map[detname+"_Y"] = readoutboardSize_Y;
  dimension_map[detname+"_Z"] = readoutboardSize_Z;

  // Connectors/flex/stiffener
  detname = detbase + "_Connectors";
  double connectorSize_X = servicehybridSize_X;
  double connectorSize_Y = servicehybridSize_Y;
  double connectorSize_Z = 1.5*mm;
  dimension_map[detname+"_X"] = connectorSize_X;
  dimension_map[detname+"_Y"] = connectorSize_Y;
  dimension_map[detname+"_Z"] = connectorSize_Z;

  // Power board
  detname = detbase + "_PowerBoard";
  double powerboardSize_X = servicehybridSize_X;
  double powerboardSize_Y = servicehybridSize_Y;
  double powerboardSize_Z = 3.1*mm;
  dimension_map[detname+"_X"] = powerboardSize_X;
  dimension_map[detname+"_Y"] = powerboardSize_Y;
  dimension_map[detname+"_Z"] = powerboardSize_Z;

  detname = detbase;
  double servicehybridSize_Z = (thermalpadSize_Z+readoutboardSize_Z+connectorSize_Z+powerboardSize_Z);
  dimension_map[detname+"_X"] = servicehybridSize_X;
  dimension_map[detname+"_Y"] = servicehybridSize_Y;
  dimension_map[detname+"_Z"] = servicehybridSize_Z;
}

