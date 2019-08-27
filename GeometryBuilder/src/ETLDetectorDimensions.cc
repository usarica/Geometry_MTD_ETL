#include <cassert>
#include <map>
#include <algorithm>
#include <utility>
#include "ETLDetectorDimensions.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"


// Could use unordered_map as well
typedef std::map<std::string, double> dim_map_t;
typedef std::map<std::string, double>::iterator dim_map_it_t;
typedef std::map<std::string, double>::const_iterator dim_map_cit_t;


namespace ETLDetectorDimensions{
  dim_map_t dimension_map;

  void setETLOffsets();
  void setETLOneSensorModuleDimensions();
  void setETLTwoSensorModuleDimensions();
  void setETLSensorServiceHybridDimensions(int const& nSensorsPerSide); // nSensorsPerSide = 3 or 6
  void setETLWedgeDimensions();
  void setETLSupports();

  void printDimensions();
}


using namespace std;
using namespace CLHEP;


double ETLDetectorDimensions::getDimension(std::string const& subdet){
  dim_map_cit_t it = dimension_map.find(subdet);
  if (it != dimension_map.cend()) return it->second;
  else{
    G4cerr << "ETLDetectorDimensions::getDimension: Could not find the dimension " << subdet << G4endl;
    assert(0);
    return -1;
  }
}
void ETLDetectorDimensions::setDimensions(){
  setETLOffsets();
  setETLTwoSensorModuleDimensions();
  setETLOneSensorModuleDimensions();
  setETLSensorServiceHybridDimensions(3);
  setETLSensorServiceHybridDimensions(6);
  setETLWedgeDimensions();
  setETLSupports();

  printDimensions();
}

void ETLDetectorDimensions::printDimensions(){ for (auto it:dimension_map) G4cout << "Dimension " << it.first << " = " << it.second << G4endl; }

void ETLDetectorDimensions::setETLOneSensorModuleDimensions(){
  string const detbase = "ETLOneSensorModule";
  string detname;

  // AlN base plate
  detname = detbase + "_BasePlate";
  const double baseplateSize_X = 28.25*mm;
  const double baseplateSize_Y = 43.1*mm;
  const double baseplateSize_Z = 0.79*mm;
  dimension_map[detname+"_X"] = baseplateSize_X;
  dimension_map[detname+"_Y"] = baseplateSize_Y;
  dimension_map[detname+"_Z"] = baseplateSize_Z;


  // Thermal pad
  detname = detbase + "_BaseFilm";
  const double basefilmSize_X = baseplateSize_X;
  const double basefilmSize_Y = baseplateSize_Y;
  const double basefilmSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = basefilmSize_X;
  dimension_map[detname+"_Y"] = basefilmSize_Y;
  dimension_map[detname+"_Z"] = basefilmSize_Z;

  // ETROC
  detname = detbase + "_ETROC";
  const double etrocSize_X = 22.3*mm;
  const double etrocSize_Y = 20.8*mm;
  const double etrocSize_Z = 0.25*mm;
  const double etrocOffset_X = getDimension("ETLTwoSensorModule_ETROC_Sep_X")/2.;
  const double etrocSep_Y = getDimension("ETLTwoSensorModule_ETROC_Sep_Y");
  dimension_map[detname+"_X"] = etrocSize_X;
  dimension_map[detname+"_Y"] = etrocSize_Y;
  dimension_map[detname+"_Z"] = etrocSize_Z;
  dimension_map[detname+"_Offset_X"] = etrocOffset_X;
  dimension_map[detname+"_Sep_Y"] = etrocSep_Y;
  if (baseplateSize_X<etrocSize_X+etrocOffset_X){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: Base plate X size < ETROC X size + ETROC X offset!" << G4endl;
    assert(0);
  }
  if (baseplateSize_Y<etrocSize_Y*2.+etrocSep_Y){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: Base plate Y size < ETROC Y size x 2 + ETROC Y separation!" << G4endl;
    assert(0);
  }

  // Laird film
  detname = detbase + "_LairdFilm";
  const double lairdfilmSize_X = etrocSize_X;
  const double lairdfilmSize_Y = etrocSize_Y*2.+etrocSep_Y;
  const double lairdfilmSize_Z = 0.08*mm;
  const double lairdfilmOffset_X = etrocOffset_X;
  dimension_map[detname+"_X"] = lairdfilmSize_X;
  dimension_map[detname+"_Y"] = lairdfilmSize_Y;
  dimension_map[detname+"_Z"] = lairdfilmSize_Z;
  dimension_map[detname+"_Offset_X"] = lairdfilmOffset_X;

  // LGAD sensor
  detname = detbase + "_LGAD";
  const double lgadSize_X = 21.2*mm;
  const double lgadSize_Y = 42.0*mm;
  const double lgadSize_Z = 0.3*mm;
  const double lgadOffset_X = getDimension("ETLTwoSensorModule_LGAD_Sep_X")/2.;
  dimension_map[detname+"_X"] = lgadSize_X;
  dimension_map[detname+"_Y"] = lgadSize_Y;
  dimension_map[detname+"_Z"] = lgadSize_Z;
  dimension_map[detname+"_Offset_X"] = lgadOffset_X;
  if (lgadSize_Y<etrocSize_Y*2.+etrocSep_Y){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: LGAD Y size < ETROC Y size x 2 + ETROC Y separation!" << G4endl;
    assert(0);
  }
  if (baseplateSize_X<lgadSize_X+lgadOffset_X){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: Base plate X size < LGAD X size + LGAD X offset!" << G4endl;
    assert(0);
  }

  // Solder bumps
  detname = detbase + "_SolderBumps";
  const double bumpsSize_X = lgadSize_X;
  const double bumpsSize_Y = etrocSize_Y;
  const double bumpsSize_Z = 0.03*mm;
  const double bumpsOffset_X = lgadOffset_X;
  const double bumpsSep_Y = etrocSep_Y;
  dimension_map[detname+"_X"] = bumpsSize_X;
  dimension_map[detname+"_Y"] = bumpsSize_Y;
  dimension_map[detname+"_Z"] = bumpsSize_Z;
  dimension_map[detname+"_Offset_X"] = bumpsOffset_X;
  dimension_map[detname+"_Sep_Y"] = bumpsSep_Y;

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  const double epoxySize_X = lgadSize_X;
  const double epoxySize_Y = lgadSize_Y;
  const double epoxySize_Z = 0.08*mm;
  const double epoxyOffset_X = lgadOffset_X;
  dimension_map[detname+"_X"] = epoxySize_X;
  dimension_map[detname+"_Y"] = epoxySize_Y;
  dimension_map[detname+"_Z"] = epoxySize_Z;
  dimension_map[detname+"_Offset_X"] = epoxyOffset_X;

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  const double coverplateSize_X = lgadSize_X; // Slightly wrong, should be slightly longer
  const double coverplateSize_Y = lgadSize_Y;
  const double coverplateSize_Z = 0.51*mm;
  const double coverplateOffset_X = getDimension("ETLTwoSensorModule_CoverPlate_Sep_X")/2.;
  dimension_map[detname+"_X"] = coverplateSize_X;
  dimension_map[detname+"_Y"] = coverplateSize_Y;
  dimension_map[detname+"_Z"] = coverplateSize_Z;
  dimension_map[detname+"_Offset_X"] = coverplateOffset_X;

  // All of the module
  detname = detbase;
  const double moduleSize_X = baseplateSize_X;
  const double moduleSize_Y = baseplateSize_Y;
  const double moduleSize_Z = basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+epoxySize_Z+lgadSize_Z+coverplateSize_Z;
  dimension_map[detname+"_X"] = moduleSize_X;
  dimension_map[detname+"_Y"] = moduleSize_Y;
  dimension_map[detname+"_Z"] = moduleSize_Z;
}
void ETLDetectorDimensions::setETLTwoSensorModuleDimensions(){
  string const detbase = "ETLTwoSensorModule";
  string detname;

  // AlN base plate
  detname = detbase + "_BasePlate";
  const double baseplateSize_X = 56.5*mm;
  const double baseplateSize_Y = 43.1*mm;
  const double baseplateSize_Z = 0.79*mm;
  dimension_map[detname+"_X"] = baseplateSize_X;
  dimension_map[detname+"_Y"] = baseplateSize_Y;
  dimension_map[detname+"_Z"] = baseplateSize_Z;

  // Thermal pad
  detname = detbase + "_BaseFilm";
  const double basefilmSize_X = baseplateSize_X;
  const double basefilmSize_Y = baseplateSize_Y;
  const double basefilmSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = basefilmSize_X;
  dimension_map[detname+"_Y"] = basefilmSize_Y;
  dimension_map[detname+"_Z"] = basefilmSize_Z;

  // ETROC
  detname = detbase + "_ETROC";
  const double etrocSize_X = 22.3*mm;
  const double etrocSize_Y = 20.8*mm;
  const double etrocSize_Z = 0.25*mm;
  const double etrocSep_X = 0.3*mm;
  const double etrocSep_Y = 0.1*mm;
  dimension_map[detname+"_X"] = etrocSize_X;
  dimension_map[detname+"_Y"] = etrocSize_Y;
  dimension_map[detname+"_Z"] = etrocSize_Z;
  dimension_map[detname+"_Sep_X"] = etrocSep_X;
  dimension_map[detname+"_Sep_Y"] = etrocSep_Y;
  if (baseplateSize_X<etrocSize_X*2.+etrocSep_X){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: Base plate X size < ETROC X size x 2 + ETROC X separation!" << G4endl;
    assert(0);
  }
  if (baseplateSize_Y<etrocSize_Y*2.+etrocSep_Y){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: Base plate Y size < ETROC Y size x 2 + ETROC Y separation!" << G4endl;
    assert(0);
  }

  // Laird film
  detname = detbase + "_LairdFilm";
  const double lairdfilmSize_X = etrocSize_X*2.+etrocSep_X;
  const double lairdfilmSize_Y = etrocSize_Y*2.+etrocSep_Y;
  const double lairdfilmSize_Z = 0.08*mm;
  dimension_map[detname+"_X"] = lairdfilmSize_X;
  dimension_map[detname+"_Y"] = lairdfilmSize_Y;
  dimension_map[detname+"_Z"] = lairdfilmSize_Z;

  // LGAD sensor
  detname = detbase + "_LGAD";
  const double lgadSize_X = 21.2*mm;
  const double lgadSize_Y = 42.0*mm;
  const double lgadSize_Z = 0.3*mm;
  const double lgadSep_X = 0.25*mm;
  dimension_map[detname+"_X"] = lgadSize_X;
  dimension_map[detname+"_Y"] = lgadSize_Y;
  dimension_map[detname+"_Z"] = lgadSize_Z;
  dimension_map[detname+"_Sep_X"] = lgadSep_X;
  if (lgadSize_Y<etrocSize_Y*2.+etrocSep_Y){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: LGAD Y size < ETROC Y size x 2 + ETROC Y separation!" << G4endl;
    assert(0);
  }
  if (baseplateSize_X<lgadSize_X*2.+lgadSep_X){
    G4cerr << "ETLDetectorDimensions::setETLTwoSensorModuleDimensions: Base plate X size < LGAD X size x 2 + LGAD X separation!" << G4endl;
    assert(0);
  }

  // Solder bumps
  detname = detbase + "_SolderBumps";
  const double bumpsSize_X = lgadSize_X;
  const double bumpsSize_Y = etrocSize_Y;
  const double bumpsSize_Z = 0.03*mm;
  const double bumpsSep_X = lgadSep_X;
  const double bumpsSep_Y = etrocSep_Y;
  dimension_map[detname+"_X"] = bumpsSize_X;
  dimension_map[detname+"_Y"] = bumpsSize_Y;
  dimension_map[detname+"_Z"] = bumpsSize_Z;
  dimension_map[detname+"_Sep_X"] = bumpsSep_X;
  dimension_map[detname+"_Sep_Y"] = bumpsSep_Y;

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  const double epoxySize_X = lgadSize_X;
  const double epoxySize_Y = lgadSize_Y;
  const double epoxySize_Z = 0.08*mm;
  const double epoxySep_X = lgadSep_X;
  dimension_map[detname+"_X"] = epoxySize_X;
  dimension_map[detname+"_Y"] = epoxySize_Y;
  dimension_map[detname+"_Z"] = epoxySize_Z;
  dimension_map[detname+"_Sep_X"] = epoxySep_X;

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  const double coverplateSize_X = lgadSize_X; // FIXME: Slightly wrong, should be slightly longer
  const double coverplateSize_Y = lgadSize_Y;
  const double coverplateSize_Z = 0.51*mm;
  const double coverplateSep_X = lgadSep_X; // FIXME: Separation is different from LGADs based on Fig. 3.73
  dimension_map[detname+"_X"] = coverplateSize_X;
  dimension_map[detname+"_Y"] = coverplateSize_Y;
  dimension_map[detname+"_Z"] = coverplateSize_Z;
  dimension_map[detname+"_Sep_X"] = coverplateSep_X;

  // All of the module
  detname = detbase;
  const double moduleSize_X = baseplateSize_X;
  const double moduleSize_Y = baseplateSize_Y;
  const double moduleSize_Z = basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+epoxySize_Z+lgadSize_Z+coverplateSize_Z;
  dimension_map[detname+"_X"] = moduleSize_X;
  dimension_map[detname+"_Y"] = moduleSize_Y;
  dimension_map[detname+"_Z"] = moduleSize_Z;
}
void ETLDetectorDimensions::setETLSensorServiceHybridDimensions(int const& nSensorsPerSide){
  string const detbase = "ETL" + std::to_string(2*nSensorsPerSide) + "SensorServiceHybrid";
  string detname;

  const double servicehybridSize_X = 34*mm;
  const double servicehybridSize_Y = getDimension("ETLTwoSensorModule_BasePlate_Y")*((double) nSensorsPerSide) + ((double) (nSensorsPerSide-1))*getDimension("ETLOffset_Module_Module_dY");

  // Thermal pad
  detname = detbase + "_ThermalPad";
  const double thermalpadSize_X = servicehybridSize_X;
  const double thermalpadSize_Y = servicehybridSize_Y;
  const double thermalpadSize_Z = 0.25*mm;
  dimension_map[detname+"_X"] = thermalpadSize_X;
  dimension_map[detname+"_Y"] = thermalpadSize_Y;
  dimension_map[detname+"_Z"] = thermalpadSize_Z;

  // Readout board
  detname = detbase + "_ReadoutBoard";
  const double readoutboardSize_X = servicehybridSize_X;
  const double readoutboardSize_Y = servicehybridSize_Y;
  const double readoutboardSize_Z = 1.*mm;
  dimension_map[detname+"_X"] = readoutboardSize_X;
  dimension_map[detname+"_Y"] = readoutboardSize_Y;
  dimension_map[detname+"_Z"] = readoutboardSize_Z;

  // Connectors/flex/stiffener
  detname = detbase + "_Connectors";
  const double connectorSize_X = servicehybridSize_X;
  const double connectorSize_Y = servicehybridSize_Y;
  const double connectorSize_Z = 1.5*mm;
  dimension_map[detname+"_X"] = connectorSize_X;
  dimension_map[detname+"_Y"] = connectorSize_Y;
  dimension_map[detname+"_Z"] = connectorSize_Z;

  // Power board
  detname = detbase + "_PowerBoard";
  const double powerboardSize_X = servicehybridSize_X;
  const double powerboardSize_Y = servicehybridSize_Y;
  const double powerboardSize_Z = 3.1*mm;
  dimension_map[detname+"_X"] = powerboardSize_X;
  dimension_map[detname+"_Y"] = powerboardSize_Y;
  dimension_map[detname+"_Z"] = powerboardSize_Z;

  detname = detbase;
  const double servicehybridSize_Z = (thermalpadSize_Z+readoutboardSize_Z+connectorSize_Z+powerboardSize_Z);
  dimension_map[detname+"_X"] = servicehybridSize_X;
  dimension_map[detname+"_Y"] = servicehybridSize_Y;
  dimension_map[detname+"_Z"] = servicehybridSize_Z;
}
void ETLDetectorDimensions::setETLWedgeDimensions(){
  string const detbase = "ETLWedge";
  string detname;

  // Wedge
  detname = detbase;
  const double wedge_rmin = 0.31*m;
  const double wedge_rmax = 1.27*m;
  const double wedge_MIC6Al_z = 6.35*mm;
  const double wedge_Epoxy_z = 0.08*mm;
  const double wedge_CoolingAl_z = 0.81*mm;
  const double wedge_z = wedge_MIC6Al_z + wedge_Epoxy_z + wedge_CoolingAl_z; // = 7.24 mm
  const double wedge_fullz = wedge_z + 2.*std::max(
    std::max(getDimension("ETLOneSensorModule_Z"), getDimension("ETLTwoSensorModule_Z")),
    std::max(getDimension("ETL6SensorServiceHybrid_Z"), getDimension("ETL12SensorServiceHybrid_Z"))
  ); // = 7.24 mm + 5.85 mm * 2 = 18.94 mm
  dimension_map[detname+"_Rmin"] = wedge_rmin;
  dimension_map[detname+"_Rmax"] = wedge_rmax;
  dimension_map[detname+"_MIC6Al_Z"] = wedge_MIC6Al_z;
  dimension_map[detname+"_Epoxy_Z"] = wedge_Epoxy_z;
  dimension_map[detname+"_CoolingAl_Z"] = wedge_CoolingAl_z;
  dimension_map[detname+"_Z"] = wedge_z; // This is the z extent of the wedge components.
  dimension_map[detname+"_FullZ"] = wedge_fullz; // This is the full z extent of the wedge including modules and service hybrids.

  detname = detbase + "_CoolingPipe";
  const double coolingpipe_rmin = 4.*mm/2.; // FIXME: Not sure what the actual value should (will?) be
  const double coolingpipe_rmax = 4.2*mm/2.; // FIXME: Not sure what the actual value should (will?) be
  dimension_map[detname+"_Rmin"] = coolingpipe_rmin;
  dimension_map[detname+"_Rmax"] = coolingpipe_rmax;
  if (coolingpipe_rmin>=coolingpipe_rmax){
    G4cerr << "ETLDetectorDimensions::setETLWedgeDimensions: The cooling pipe has to have an outer radius greater than its inner radius!" << G4endl;
    assert(0);
  }
  if (coolingpipe_rmax*2.>wedge_MIC6Al_z){
    G4cerr << "ETLDetectorDimensions::setETLWedgeDimensions: The outer diameter of the cooling pipe is larger than the thickness of the MIC6-aluminum wedge component!" << G4endl;
    assert(0);
  }
  if (coolingpipe_rmax>=wedge_rmax){
    G4cerr << "ETLDetectorDimensions::setETLWedgeDimensions: The outer radius of the cooling pipe has to be smaller than the outer radius of the wedge!" << G4endl;
    assert(0);
  }

  detname = detbase + "_Attachment";
  const double attachment_x = (wedge_rmax - wedge_rmin);
  const double attachment_y = 5.*cm;
  const double attachment_z = 1.*mm;
  const double attachment_offset_x = wedge_rmin;
  const double attachment_offset_y = attachment_y*(0.5 - 0.5);
  dimension_map[detname+"_X"] = attachment_x;
  dimension_map[detname+"_Y"] = attachment_y;
  dimension_map[detname+"_Z"] = attachment_z;
  dimension_map[detname+"_Offset_X"] = attachment_offset_x;
  dimension_map[detname+"_Offset_Y"] = attachment_offset_y;
  if (attachment_z>(wedge_fullz-wedge_z)/2.){
    G4cerr << "ETLDetectorDimensions::setETLWedgeDimensions: Wedge front attachment support bar has a dz value greater than what is allowed by the full dz of the wedge!" << G4endl;
    assert(0);
  }
  if (attachment_y>attachment_x){
    G4cerr << "ETLDetectorDimensions::setETLWedgeDimensions: Wedge front attachment support bar has a dy value greater than its dx (= wedge delta R) value!" << G4endl;
    assert(0);
  }

}
void ETLDetectorDimensions::setETLOffsets(){
  // This function only handles inter-component spacing
  // Offsets within two-sensor modules are supposed to be inside setETLTwoSensorModuleDimensions.

  string const detbase = "ETLOffset";
  string detname;

  // dX between any module and a service hybrid
  detname = detbase + "_Module_SensorServiceHybrid_dX";
  const double sep_X_module_servicehybrid = 0.*mm;
  dimension_map[detname] = sep_X_module_servicehybrid;

  // dY between any module
  detname = detbase + "_Module_Module_dY";
  const double sep_Y_module_module = 0.5*mm;
  dimension_map[detname] = sep_Y_module_module;

  // dX between first module and wedge attachment
  detname = detbase + "_Module_WedgeAttachment_dX";
  const double sep_X_module_wedgeAttachment = sep_Y_module_module;
  dimension_map[detname] = sep_X_module_wedgeAttachment;

  // dY between first module and wedge attachment
  detname = detbase + "_Module_WedgeAttachment_dY";
  const double sep_Y_module_wedgeAttachment = sep_Y_module_module;
  dimension_map[detname] = sep_Y_module_wedgeAttachment;

  // dY between any service hybrid
  detname = detbase + "_SensorServiceHybrid_SensorServiceHybrid_dY";
  const double sep_Y_servicehybrid_servicehybrid = sep_Y_module_module;
  dimension_map[detname] = sep_Y_servicehybrid_servicehybrid;

  // dZ between two disks
  detname = detbase + "_Disks_dZ";
  const double sep_Z_disks = 20.*mm;
  dimension_map[detname] = sep_Z_disks;

  // ~dZ of ETL CoM
  detname = detbase + "_IP_dZ";
  const double IP_dZ = 2.98*m;
  dimension_map[detname] = IP_dZ;
}
void ETLDetectorDimensions::setETLSupports(){
  // This function handles the support structure dimensions.

  string const detbase = "ETLSupport";
  string detname;

  // Polyethlene neutron moderator
  detname = detbase + "_NeutronModerator";
  const double neutronModerator_Z = 12.*cm;
  dimension_map[detname+"_Z"] = neutronModerator_Z;

  // Support tube
  detname = detbase + "_SupportTube_Thickness";
  const double supportTube_thickness = 6.*mm;
  dimension_map[detname] = supportTube_thickness;

  // Mounting bracket
  detname = detbase + "_MountingBracket_Thickness";
  const double mountingBracket_thickness = supportTube_thickness;
  dimension_map[detname] = mountingBracket_thickness;

  // Thermal screen (aerogel core) permaglas skin thickness
  detname = detbase + "_ThermalScreen_Skin_Thickness";
  const double thermalScreen_skin_thickness = 1.5*mm;
  dimension_map[detname] = thermalScreen_skin_thickness;
}
