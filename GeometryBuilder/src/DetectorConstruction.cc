#include <cassert>
#include <unordered_map>
#include "DetectorConstruction.hh"
#include "ETLDetectorDimensions.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4ReflectionFactory.hh"
#include "G4SystemOfUnits.hh"


namespace DetectorConstructionHelpers{
  BasicDetectorAttributes::BasicDetectorAttributes() :
    solid(nullptr),
    logical(nullptr),
    material(nullptr),
    visualization(nullptr)
  {}
  BasicDetectorAttributes::BasicDetectorAttributes(BasicDetectorAttributes const& other) :
    solid(other.solid),
    logical(other.logical),
    material(other.material),
    visualization(other.visualization)
  {}
  BasicDetectorAttributes::BasicDetectorAttributes(G4VSolid* solid_, G4LogicalVolume* logical_, G4Material* material_, G4VisAttributes* visualization_) :
    solid(solid_),
    logical(logical_),
    material(material_),
    visualization(visualization_)
  {}

  std::unordered_map<std::string, BasicDetectorAttributes> detector_components;
}


using namespace CLHEP;
using namespace std;
using namespace DetectorConstructionHelpers;


DetectorConstruction::DetectorConstruction() :
  G4VUserDetectorConstruction(),
  fCheckOverlaps(true),
  fScoringVolume(nullptr)
{}
DetectorConstruction::~DetectorConstruction(){}


G4VPhysicalVolume* DetectorConstruction::Construct(){
  G4cout << "Begin DetectorConstruction::Construct" << G4endl;
  DefineMaterials();
  ETLDetectorDimensions::setDimensions(); // Build the dimensions
  G4VPhysicalVolume* res = DefineVolumes();
  G4cout << "End DetectorConstruction::Construct" << G4endl;
  return res;
}
void DetectorConstruction::DefineMaterials(){
  G4cout << "Begin DetectorConstruction::DefineMaterials" << G4endl;

  G4NistManager* nist = G4NistManager::Instance();
  nist->ListMaterials("all");
  G4double a;
  G4double z;
  G4double density;

  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* nist_CO2 = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4Material* nist_Al = nist->FindOrBuildMaterial("G4_Al");
  G4Material* nist_Sn = nist->FindOrBuildMaterial("G4_Sn");
  G4Material* nist_Ag = nist->FindOrBuildMaterial("G4_Ag");
  G4Material* nist_Cu = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* nist_Steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material* nist_SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* nist_Glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

  nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  //G4double env_temp = (273.16+15.00)*kelvin; // NTP_Temperature=293.15*kelvin
  //G4double env_pressure = STP_Pressure; // in units of pascal

  // Vacuum
  new G4Material(
    "Vacuum", z=1., a=1.01*g/mole, density=universe_mean_density,
    kStateGas, 2.73*kelvin, 3.e-18*pascal
  );


  G4cout << "DetectorConstruction::DefineMaterials: Environment materials built." << G4endl;

  // MIC6 Aluminum
  // Used as ETL disk material
  nist->BuildMaterialWithNewDensity("MIC6_Al", "G4_Al", 2.796*g/cm3, NTP_Temperature, STP_Pressure); // 0.101*lb/in3
  nist->BuildMaterialWithNewDensity("Cool_CO2", "G4_CARBON_DIOXIDE", nist_CO2->GetDensity() * NTP_Temperature / ((273.15 - 20.-35.)*CLHEP::kelvin), (273.15 - 20.-35.)*CLHEP::kelvin, STP_Pressure);


  // When calling FindOrBuildElement, do not use 'G4_' at the beginning of the element name
  G4Element* Aluminum = nist->FindOrBuildElement("Al");
  G4Element* Nitrogen = nist->FindOrBuildElement("N");
  G4Element* Carbon = nist->FindOrBuildElement("C");
  G4Element* Hydrogen = nist->FindOrBuildElement("H");
  G4Element* Oxygen = nist->FindOrBuildElement("O");
  G4Element* Silicon = nist->FindOrBuildElement("Si");
  G4Element* Magnesium = nist->FindOrBuildElement("Mg");
  //G4Element* Tin = nist->FindOrBuildElement("Sn");
  //G4Element* Silver = nist->FindOrBuildElement("Ag");
  G4Element* Copper = nist->FindOrBuildElement("Cu");
  G4cout << "DetectorConstruction::DefineMaterials: Elements built." << G4endl;

  if (!Aluminum) G4cerr << "DetectorConstruction::DefineMaterials: Element Al not found!" << G4endl;
  if (!Nitrogen) G4cerr << "DetectorConstruction::DefineMaterials: Element N not found!" << G4endl;
  if (!Carbon) G4cerr << "DetectorConstruction::DefineMaterials: Element C not found!" << G4endl;
  if (!Silicon) G4cerr << "DetectorConstruction::DefineMaterials: Element Si not found!" << G4endl;

  // Aluminum Nitride (AlN), used in sensors and elsewhere
  G4Material* AlN = new G4Material("AlN", 3.260*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure);
  AlN->AddElement(Aluminum, (G4int) 1);
  AlN->AddElement(Nitrogen, (G4int) 1);

  G4Material* Mg_OH2 = new G4Material("Mg_OH2", 2.34*g/cm3, 3, kStateSolid, NTP_Temperature, STP_Pressure);
  Mg_OH2->AddElement(Magnesium, (G4int) 1);
  Mg_OH2->AddElement(Hydrogen, (G4int) 2);
  Mg_OH2->AddElement(Oxygen, (G4int) 2);

  G4Material* PET = new G4Material("PET", 1.38*g/cm3, 3, kStateSolid, NTP_Temperature, STP_Pressure);
  PET->AddElement(Carbon, (G4int) 10);
  PET->AddElement(Hydrogen, (G4int) 8);
  PET->AddElement(Oxygen, (G4int) 4);

  G4Material* FR4 = new G4Material("FR4", 1.85*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure); // For circuit boards, density approximate?
  FR4->AddElement(Carbon, (G4int) 1);
  FR4->AddElement(Hydrogen, (G4int) 1);

  G4Material* PCB = new G4Material("PCB", 1.85*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure); // For circuit boards, density approximate?
  PCB->AddElement(Carbon, (G4int) 1);
  PCB->AddElement(Hydrogen, (G4int) 1);

  G4Material* Permaglas = new G4Material("Permaglas", 1.85*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure); // Insulation material, density approximate
  Permaglas->AddElement(Carbon, (G4int) 1);
  Permaglas->AddElement(Hydrogen, (G4int) 1);

  G4Material* Mat_LGAD = new G4Material("Mat_LGAD", 2.3290*g/cm3, 1, kStateSolid, NTP_Temperature, STP_Pressure); // For LGAD or ETROC
  Mat_LGAD->AddElement(Silicon, (G4int) 1);

  G4Material* Mat_ETROC = new G4Material("Mat_ETROC", 2.3290*g/cm3, 1, kStateSolid, NTP_Temperature, STP_Pressure); // For LGAD or ETROC
  Mat_ETROC->AddElement(Silicon, (G4int) 1);

  G4Material* Epoxy = new G4Material("Epoxy", 1.1*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure); // Density approximate
  Epoxy->AddElement(Carbon, (G4int) 1);
  Epoxy->AddElement(Hydrogen, (G4int) 1);

  G4Material* Laird = new G4Material("Laird", 1.*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure); // Density approximate
  Laird->AddElement(Carbon, (G4int) 1);
  Laird->AddElement(Hydrogen, (G4int) 1);

  G4Material* SnAg = new G4Material("SnAg", (nist_Sn->GetDensity()*0.5 + nist_Ag->GetDensity()*0.5)*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure); // FIXME: Composition should be per mol, and density needs recalc.
  SnAg->AddMaterial(nist_Sn, (G4double) 0.5);
  SnAg->AddMaterial(nist_Ag, (G4double) 0.5);

  G4Material* SensorBump = new G4Material("SensorBump", (SnAg->GetDensity()*0.1 + Air->GetDensity()*0.9), 2, kStateSolid, NTP_Temperature, STP_Pressure); // FIXME: Composition should be per mol, and density needs recalc
  SensorBump->AddMaterial(SnAg, (G4double) 0.1);
  SensorBump->AddMaterial(Air, (G4double) 0.9);

  G4Material* ServiceConnector = new G4Material("ServiceConnector", (nist_Cu->GetDensity()*0.1 + Air->GetDensity()*0.9), 2, kStateSolid, NTP_Temperature, STP_Pressure); // FIXME: Composition should be per mol, and density needs recalc
  ServiceConnector->AddMaterial(nist_Cu, (G4double) 0.1);
  ServiceConnector->AddMaterial(Air, (G4double) 0.9);

  G4Material* Aerogel = new G4Material("Aerogel", 0.16*g/cm3, 6, kStateSolid, NTP_Temperature, STP_Pressure);
  Aerogel->AddMaterial(Air, (G4double) 0.9);
  Aerogel->AddMaterial(nist_SiO2, (G4double) 0.065);
  Aerogel->AddMaterial(PET, (G4double) 0.015);
  Aerogel->AddMaterial(nist_Glass, (G4double) 0.015);
  Aerogel->AddMaterial(Mg_OH2, (G4double) 0.0025);
  Aerogel->AddMaterial(nist_Al, (G4double) 0.0025);

  G4cout << "DetectorConstruction::DefineMaterials: AlN built!" << G4endl;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


void DetectorConstruction::BuildOneSensorModule(
  G4LogicalVolume* motherLogical,
  G4RotationMatrix* rotation, G4ThreeVector const& relativePos
){
  using namespace ETLDetectorDimensions;

  // See Fig. 3.56 in the MTD TDR for the module diagrams
  // Also Fig. 3.73 for the z positions.
  // Note that the base aluminum plate on which the service is loaded, or the epoxy underneath are not included since they are only on one side.
  string const detbase = "ETLOneSensorModule";
  string detname = detbase;

  // Search for the module
  auto det_component = detector_components.find(detname);
  if (det_component!=detector_components.end()){
    new G4PVPlacement(
      rotation, relativePos,
      det_component->second.logical,
      detname.c_str(),
      motherLogical,
      false, 0, fCheckOverlaps
    );
    return;
  }

  // AlN base plate
  detname = detbase + "_BasePlate";
  G4double baseplateSize_X = getDimension(detname+"_X");
  G4double baseplateSize_Y = getDimension(detname+"_Y");
  G4double baseplateSize_Z = getDimension(detname+"_Z");
  G4Material* baseplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* baseplateVisAttr = new G4VisAttributes(G4Colour::Gray()); baseplateVisAttr->SetVisibility(true);

  // Thermal pad
  detname = detbase + "_BaseFilm";
  G4double basefilmSize_X = getDimension(detname+"_X");
  G4double basefilmSize_Y = getDimension(detname+"_Y");
  G4double basefilmSize_Z = getDimension(detname+"_Z");
  G4Material* basefilmMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* basefilmVisAttr = new G4VisAttributes(G4Colour::Brown()); basefilmVisAttr->SetVisibility(true);

  // ETROC
  detname = detbase + "_ETROC";
  G4double etrocSize_X = getDimension(detname+"_X");
  G4double etrocSize_Y = getDimension(detname+"_Y");
  G4double etrocSize_Z = getDimension(detname+"_Z");
  G4double etrocOffset_X = getDimension(detname+"_Offset_X");
  G4double etrocSep_Y = getDimension(detname+"_Sep_Y");
  G4Material* etrocMat = G4Material::GetMaterial("Mat_ETROC");
  G4VisAttributes* etrocVisAttr = new G4VisAttributes(G4Colour::Green()); etrocVisAttr->SetVisibility(true);

  // Laird film
  detname = detbase + "_LairdFilm";
  G4double lairdfilmSize_X = getDimension(detname+"_X");
  G4double lairdfilmSize_Y = getDimension(detname+"_Y");
  G4double lairdfilmSize_Z = getDimension(detname+"_Z");
  G4double lairdfilmOffset_X = getDimension(detname+"_Offset_X");
  G4Material* lairdfilmMat = G4Material::GetMaterial("Laird");
  G4VisAttributes* lairdfilmVisAttr = new G4VisAttributes(G4Colour::Brown()); lairdfilmVisAttr->SetVisibility(true);

  // LGAD sensor
  detname = detbase + "_LGAD";
  G4double lgadSize_X = getDimension(detname+"_X");
  G4double lgadSize_Y = getDimension(detname+"_Y");
  G4double lgadSize_Z = getDimension(detname+"_Z");
  G4double lgadOffset_X = getDimension(detname+"_Offset_X");
  G4Material* lgadMat = G4Material::GetMaterial("Mat_LGAD");
  G4VisAttributes* lgadVisAttr = new G4VisAttributes(G4Colour::Yellow()); lgadVisAttr->SetVisibility(true);

  // Solder bumps
  detname = detbase + "_SolderBumps";
  G4double bumpsSize_X = getDimension(detname+"_X");
  G4double bumpsSize_Y = getDimension(detname+"_Y");
  G4double bumpsSize_Z = getDimension(detname+"_Z");
  G4double bumpsOffset_X = getDimension(detname+"_Offset_X");
  G4double bumpsSep_Y = getDimension(detname+"_Sep_Y");
  G4Material* bumpsMat = G4Material::GetMaterial("SensorBump");
  G4VisAttributes* bumpsVisAttr = new G4VisAttributes(G4Colour::Magenta()); bumpsVisAttr->SetVisibility(true);

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  G4double epoxySize_X = getDimension(detname+"_X");
  G4double epoxySize_Y = getDimension(detname+"_Y");
  G4double epoxySize_Z = getDimension(detname+"_Z");
  G4double epoxyOffset_X = getDimension(detname+"_Offset_X");
  G4Material* epoxyMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* epoxyVisAttr = new G4VisAttributes(G4Colour::Brown()); epoxyVisAttr->SetVisibility(true);

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  G4double coverplateSize_X = getDimension(detname+"_X");
  G4double coverplateSize_Y = getDimension(detname+"_Y");
  G4double coverplateSize_Z = getDimension(detname+"_Z");
  G4double coverplateOffset_X = getDimension(detname+"_Offset_X");
  G4Material* coverplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* coverplateVisAttr = new G4VisAttributes(G4Colour::Gray()); coverplateVisAttr->SetVisibility(true);

  detname = detbase;
  G4double moduleSize_X = getDimension(detname+"_X");
  G4double moduleSize_Y = getDimension(detname+"_Y");
  G4double moduleSize_Z = getDimension(detname+"_Z");
  G4Material* moduleMat = G4Material::GetMaterial("Vacuum");
  G4VisAttributes* moduleVisAttr = new G4VisAttributes(G4Colour::Black()); moduleVisAttr->SetVisibility(false);

  // Build the module
  detname = detbase;
  G4Box* moduleBox = new G4Box(
    detname.c_str(),
    moduleSize_X/2., moduleSize_Y/2., moduleSize_Z/2.
  );
  G4LogicalVolume* moduleLogical = new G4LogicalVolume(
    moduleBox, 
    moduleMat,
    detname.c_str()
  );
  moduleLogical->SetVisAttributes(moduleVisAttr);
  new G4PVPlacement(
    rotation, relativePos,
    moduleLogical,
    detname.c_str(),
    motherLogical,
    false, 0, fCheckOverlaps
  );
  detector_components[detname] = BasicDetectorAttributes(moduleBox, moduleLogical, moduleMat, moduleVisAttr);

  // Build the base film at the very bottom of the module
  detname = detbase + "_BaseFilm";
  G4Box* basefilmBox = new G4Box(
    detname.c_str(),
    basefilmSize_X/2., basefilmSize_Y/2., basefilmSize_Z/2.
  );
  G4LogicalVolume* basefilmLogical = new G4LogicalVolume(
    basefilmBox,
    basefilmMat,
    detname.c_str()
  );
  basefilmLogical->SetVisAttributes(basefilmVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, (basefilmSize_Z-moduleSize_Z)/2.),
    basefilmLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the base Al plate on top of the film
  detname = detbase + "_BasePlate";
  G4Box* baseplateBox = new G4Box(
    detname.c_str(),
    baseplateSize_X/2., baseplateSize_Y/2., baseplateSize_Z/2.
  );
  G4LogicalVolume* baseplateLogical = new G4LogicalVolume(
    baseplateBox,
    baseplateMat,
    detname.c_str()
  );
  baseplateLogical->SetVisAttributes(baseplateVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, basefilmSize_Z+(baseplateSize_Z-moduleSize_Z)/2.),
    baseplateLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the laird film on top of the Al plate
  detname = detbase + "_LairdFilm";
  G4Box* lairdfilmBox = new G4Box(
    detname.c_str(),
    lairdfilmSize_X/2., lairdfilmSize_Y/2., lairdfilmSize_Z/2.
  );
  G4LogicalVolume* lairdfilmLogical = new G4LogicalVolume(
    lairdfilmBox,
    lairdfilmMat,
    detname.c_str()
  );
  lairdfilmLogical->SetVisAttributes(lairdfilmVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(((lairdfilmSize_X-baseplateSize_X)/2. + lairdfilmOffset_X), 0, basefilmSize_Z+baseplateSize_Z+(lairdfilmSize_Z-moduleSize_Z)/2.),
    lairdfilmLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the ETROCs
  detname = detbase + "_ETROC";
  G4Box* etrocBox = new G4Box(
    detname.c_str(),
    etrocSize_X/2., etrocSize_Y/2., etrocSize_Z/2.
  );
  G4LogicalVolume* etrocLogical = new G4LogicalVolume(
    etrocBox,
    etrocMat,
    detname.c_str()
  );
  etrocLogical->SetVisAttributes(etrocVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(((etrocSize_X-baseplateSize_X)/2. + etrocOffset_X), -(etrocSize_Y+etrocSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(((etrocSize_X-baseplateSize_X)/2. + etrocOffset_X), +(etrocSize_Y+etrocSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );


  // Build the solder bumps
  detname = detbase + "_SolderBumps";
  G4Box* bumpsBox = new G4Box(
    detname.c_str(),
    bumpsSize_X/2., bumpsSize_Y/2., bumpsSize_Z/2.
  );
  G4LogicalVolume* bumpsLogical = new G4LogicalVolume(
    bumpsBox,
    bumpsMat,
    detname.c_str()
  );
  bumpsLogical->SetVisAttributes(bumpsVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(((bumpsSize_X-baseplateSize_X)/2. + bumpsOffset_X), -(bumpsSize_Y+bumpsSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(((bumpsSize_X-baseplateSize_X)/2. + bumpsOffset_X), +(bumpsSize_Y+bumpsSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the LGAD
  detname = detbase + "_LGAD";
  G4Box* lgadBox = new G4Box(
    detname.c_str(),
    lgadSize_X/2., lgadSize_Y/2., lgadSize_Z/2.
  );
  G4LogicalVolume* lgadLogical = new G4LogicalVolume(
    lgadBox,
    lgadMat,
    detname.c_str()
  );
  lgadLogical->SetVisAttributes(lgadVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(((lgadSize_X-baseplateSize_X)/2. + lgadOffset_X), 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
    lgadLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the epoxy cover
  detname = detbase + "_EpoxyCover";
  G4Box* epoxyBox = new G4Box(
    detname.c_str(),
    epoxySize_X/2., epoxySize_Y/2., epoxySize_Z/2.
  );
  G4LogicalVolume* epoxyLogical = new G4LogicalVolume(
    epoxyBox,
    epoxyMat,
    detname.c_str()
  );
  epoxyLogical->SetVisAttributes(epoxyVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(((epoxySize_X-baseplateSize_X)/2. + epoxyOffset_X), 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
    epoxyLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the cover plate
  detname = detbase + "_CoverPlate";
  G4Box* coverplateBox = new G4Box(
    detname.c_str(),
    coverplateSize_X/2., coverplateSize_Y/2., coverplateSize_Z/2.
  );
  G4LogicalVolume* coverplateLogical = new G4LogicalVolume(
    coverplateBox,
    coverplateMat,
    detname.c_str()
  );
  coverplateLogical->SetVisAttributes(coverplateVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(((coverplateSize_X-baseplateSize_X)/2. + coverplateOffset_X), 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
}

void DetectorConstruction::BuildTwoSensorModule(
  G4LogicalVolume* motherLogical,
  G4RotationMatrix* rotation, G4ThreeVector const& relativePos
){
  using namespace ETLDetectorDimensions;

  // See Fig. 3.56 in the MTD TDR for the module diagrams
  // Also Fig. 3.73 for the z positions.
  // Note that the base aluminum plate on which the service is loaded, or the epoxy underneath are not included since they are only on one side.
  string const detbase = "ETLTwoSensorModule";
  string detname = detbase;

  // Search for the module
  auto det_component = detector_components.find(detname);
  if (det_component!=detector_components.end()){
    new G4PVPlacement(
      rotation, relativePos,
      det_component->second.logical,
      detname.c_str(),
      motherLogical,
      false, 0, fCheckOverlaps
    );
    return;
  }

  // AlN base plate
  detname = detbase + "_BasePlate";
  G4double baseplateSize_X = getDimension(detname+"_X");
  G4double baseplateSize_Y = getDimension(detname+"_Y");
  G4double baseplateSize_Z = getDimension(detname+"_Z");
  G4Material* baseplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* baseplateVisAttr = new G4VisAttributes(G4Colour::Gray()); baseplateVisAttr->SetVisibility(true);

  // Thermal pad
  detname = detbase + "_BaseFilm";
  G4double basefilmSize_X = getDimension(detname+"_X");
  G4double basefilmSize_Y = getDimension(detname+"_Y");
  G4double basefilmSize_Z = getDimension(detname+"_Z");
  G4Material* basefilmMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* basefilmVisAttr = new G4VisAttributes(G4Colour::Brown()); basefilmVisAttr->SetVisibility(true);

  // ETROC
  detname = detbase + "_ETROC";
  G4double etrocSize_X = getDimension(detname+"_X");
  G4double etrocSize_Y = getDimension(detname+"_Y");
  G4double etrocSize_Z = getDimension(detname+"_Z");
  G4double etrocSep_X = getDimension(detname+"_Sep_X");
  G4double etrocSep_Y = getDimension(detname+"_Sep_Y");
  G4Material* etrocMat = G4Material::GetMaterial("Mat_ETROC");
  G4VisAttributes* etrocVisAttr = new G4VisAttributes(G4Colour::Green()); etrocVisAttr->SetVisibility(true);

  // Laird film
  detname = detbase + "_LairdFilm";
  G4double lairdfilmSize_X = getDimension(detname+"_X");
  G4double lairdfilmSize_Y = getDimension(detname+"_Y");
  G4double lairdfilmSize_Z = getDimension(detname+"_Z");
  G4Material* lairdfilmMat = G4Material::GetMaterial("Laird");
  G4VisAttributes* lairdfilmVisAttr = new G4VisAttributes(G4Colour::Brown()); lairdfilmVisAttr->SetVisibility(true);

  // LGAD sensor
  detname = detbase + "_LGAD";
  G4double lgadSize_X = getDimension(detname+"_X");
  G4double lgadSize_Y = getDimension(detname+"_Y");
  G4double lgadSize_Z = getDimension(detname+"_Z");
  G4double lgadSep_X = getDimension(detname+"_Sep_X");
  G4Material* lgadMat = G4Material::GetMaterial("Mat_LGAD");
  G4VisAttributes* lgadVisAttr = new G4VisAttributes(G4Colour::Yellow()); lgadVisAttr->SetVisibility(true);

  // Solder bumps
  detname = detbase + "_SolderBumps";
  G4double bumpsSize_X = getDimension(detname+"_X");
  G4double bumpsSize_Y = getDimension(detname+"_Y");
  G4double bumpsSize_Z = getDimension(detname+"_Z");
  G4double bumpsSep_X = getDimension(detname+"_Sep_X");
  G4double bumpsSep_Y = getDimension(detname+"_Sep_Y");
  G4Material* bumpsMat = G4Material::GetMaterial("SensorBump");
  G4VisAttributes* bumpsVisAttr = new G4VisAttributes(G4Colour::Magenta()); bumpsVisAttr->SetVisibility(true);

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  G4double epoxySize_X = getDimension(detname+"_X");
  G4double epoxySize_Y = getDimension(detname+"_Y");
  G4double epoxySize_Z = getDimension(detname+"_Z");
  G4double epoxySep_X = getDimension(detname+"_Sep_X");
  G4Material* epoxyMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* epoxyVisAttr = new G4VisAttributes(G4Colour::Brown()); epoxyVisAttr->SetVisibility(true);

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  G4double coverplateSize_X = getDimension(detname+"_X");
  G4double coverplateSize_Y = getDimension(detname+"_Y");
  G4double coverplateSize_Z = getDimension(detname+"_Z");
  G4double coverplateSep_X = getDimension(detname+"_Sep_X");
  G4Material* coverplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* coverplateVisAttr = new G4VisAttributes(G4Colour::Gray()); coverplateVisAttr->SetVisibility(true);

  detname = detbase;
  G4double moduleSize_X = getDimension(detname+"_X");
  G4double moduleSize_Y = getDimension(detname+"_Y");
  G4double moduleSize_Z = getDimension(detname+"_Z");
  G4Material* moduleMat = G4Material::GetMaterial("Vacuum");
  G4VisAttributes* moduleVisAttr = new G4VisAttributes(G4Colour::Black()); moduleVisAttr->SetVisibility(false);

  // Build the module
  detname = detbase;
  G4Box* moduleBox = new G4Box(
    detname.c_str(),
    moduleSize_X/2., moduleSize_Y/2., moduleSize_Z/2.
  );
  G4LogicalVolume* moduleLogical = new G4LogicalVolume(
    moduleBox,
    moduleMat,
    detname.c_str()
  );
  moduleLogical->SetVisAttributes(moduleVisAttr);
  new G4PVPlacement(
    rotation, relativePos,
    moduleLogical,
    detname.c_str(),
    motherLogical,
    false, 0, fCheckOverlaps
  );
  detector_components[detname] = BasicDetectorAttributes(moduleBox, moduleLogical, moduleMat, moduleVisAttr);


  // Build the base film at the very bottom of the module
  detname = detbase + "_BaseFilm";
  G4Box* basefilmBox = new G4Box(
    detname.c_str(),
    basefilmSize_X/2., basefilmSize_Y/2., basefilmSize_Z/2.
  );
  G4LogicalVolume* basefilmLogical = new G4LogicalVolume(
    basefilmBox,
    basefilmMat,
    detname.c_str()
  );
  basefilmLogical->SetVisAttributes(basefilmVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, (basefilmSize_Z-moduleSize_Z)/2.),
    basefilmLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the base Al plate on top of the film
  detname = detbase + "_BasePlate";
  G4Box* baseplateBox = new G4Box(
    detname.c_str(),
    baseplateSize_X/2., baseplateSize_Y/2., baseplateSize_Z/2.
  );
  G4LogicalVolume* baseplateLogical = new G4LogicalVolume(
    baseplateBox,
    baseplateMat,
    detname.c_str()
  );
  baseplateLogical->SetVisAttributes(baseplateVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, basefilmSize_Z+(baseplateSize_Z-moduleSize_Z)/2.),
    baseplateLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the laird film on top of the Al plate
  detname = detbase + "_LairdFilm";
  G4Box* lairdfilmBox = new G4Box(
    detname.c_str(),
    lairdfilmSize_X/2., lairdfilmSize_Y/2., lairdfilmSize_Z/2.
  );
  G4LogicalVolume* lairdfilmLogical = new G4LogicalVolume(
    lairdfilmBox,
    lairdfilmMat,
    detname.c_str()
  );
  lairdfilmLogical->SetVisAttributes(lairdfilmVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, basefilmSize_Z+baseplateSize_Z+(lairdfilmSize_Z-moduleSize_Z)/2.),
    lairdfilmLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the ETROCs
  detname = detbase + "_ETROC";
  G4Box* etrocBox = new G4Box(
    detname.c_str(),
    etrocSize_X/2., etrocSize_Y/2., etrocSize_Z/2.
  );
  G4LogicalVolume* etrocLogical = new G4LogicalVolume(
    etrocBox,
    etrocMat,
    detname.c_str()
  );
  etrocLogical->SetVisAttributes(etrocVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(-(etrocSize_X+etrocSep_X)/2., -(etrocSize_Y+etrocSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(-(etrocSize_X+etrocSep_X)/2., +(etrocSize_Y+etrocSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(+(etrocSize_X+etrocSep_X)/2., -(etrocSize_Y+etrocSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"3").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(+(etrocSize_X+etrocSep_X)/2., +(etrocSize_Y+etrocSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"4").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );



  // Build the solder bumps
  detname = detbase + "_SolderBumps";
  G4Box* bumpsBox = new G4Box(
    detname.c_str(),
    bumpsSize_X/2., bumpsSize_Y/2., bumpsSize_Z/2.
  );
  G4LogicalVolume* bumpsLogical = new G4LogicalVolume(
    bumpsBox,
    bumpsMat,
    detname.c_str()
  );
  bumpsLogical->SetVisAttributes(bumpsVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(-(bumpsSize_X+bumpsSep_X)/2., -(bumpsSize_Y+bumpsSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(-(bumpsSize_X+bumpsSep_X)/2., +(bumpsSize_Y+bumpsSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(+(bumpsSize_X+bumpsSep_X)/2., -(bumpsSize_Y+bumpsSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    (detname+"3").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(+(bumpsSize_X+bumpsSep_X)/2., +(bumpsSize_Y+bumpsSep_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    (detname+"4").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the LGAD
  detname = detbase + "_LGAD";
  G4Box* lgadBox = new G4Box(
    detname.c_str(),
    lgadSize_X/2., lgadSize_Y/2., lgadSize_Z/2.
  );
  G4LogicalVolume* lgadLogical = new G4LogicalVolume(
    lgadBox,
    lgadMat,
    detname.c_str()
  );
  lgadLogical->SetVisAttributes(lgadVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(-(lgadSize_X+lgadSep_X)/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
    lgadLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(+(lgadSize_X+lgadSep_X)/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
    lgadLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the epoxy cover
  detname = detbase + "_EpoxyCover";
  G4Box* epoxyBox = new G4Box(
    detname.c_str(),
    epoxySize_X/2., epoxySize_Y/2., epoxySize_Z/2.
  );
  G4LogicalVolume* epoxyLogical = new G4LogicalVolume(
    epoxyBox,
    epoxyMat,
    detname.c_str()
  );
  epoxyLogical->SetVisAttributes(epoxyVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(-(epoxySize_X+epoxySep_X)/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
    epoxyLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(+(epoxySize_X+epoxySep_X)/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
    epoxyLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the cover plate
  detname = detbase + "_CoverPlate";
  G4Box* coverplateBox = new G4Box(
    detname.c_str(),
    coverplateSize_X/2., coverplateSize_Y/2., coverplateSize_Z/2.
  );
  G4LogicalVolume* coverplateLogical = new G4LogicalVolume(
    coverplateBox,
    coverplateMat,
    detname.c_str()
  );
  coverplateLogical->SetVisAttributes(coverplateVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(-(coverplateSize_X+coverplateSep_X)/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  new G4PVPlacement(
    0, G4ThreeVector(+(coverplateSize_X+coverplateSep_X)/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
}

void DetectorConstruction::BuildSensorServiceHybrid(
  int const& nSensorsPerSide, // 6 or 3
  G4LogicalVolume* motherLogical,
  G4RotationMatrix* rotation, G4ThreeVector const& relativePos
){
  using namespace ETLDetectorDimensions;

  string const detbase = "ETL" + std::to_string(2*nSensorsPerSide) + "SensorServiceHybrid";
  string detname = detbase;

  // Search for the module
  auto det_component = detector_components.find(detname);
  if (det_component!=detector_components.end()){
    new G4PVPlacement(
      rotation, relativePos,
      det_component->second.logical,
      detname.c_str(),
      motherLogical,
      false, 0, fCheckOverlaps
    );
    return;
  }

  // See Figs. 3.61 and 3.62 in the MTD TDR for the diagrams
  detname = detbase;
  G4double servicehybridSize_X = getDimension(detname+"_X");
  G4double servicehybridSize_Y = getDimension(detname+"_Y");
  G4double servicehybridSize_Z = getDimension(detname+"_Z");
  G4Material* servicehybridMat = G4Material::GetMaterial("Vacuum");
  G4VisAttributes* servicehybridVisAttr = new G4VisAttributes(G4Colour::Gray()); servicehybridVisAttr->SetVisibility(false);

  // Thermal pad
  detname = detbase + "_ThermalPad";
  G4double thermalpadSize_X = getDimension(detname+"_X");
  G4double thermalpadSize_Y = getDimension(detname+"_Y");
  G4double thermalpadSize_Z = getDimension(detname+"_Z");
  G4Material* thermalpadMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* thermalpadVisAttr = new G4VisAttributes(G4Colour::Brown()); thermalpadVisAttr->SetVisibility(true);

  // Readout board
  detname = detbase + "_ReadoutBoard";
  G4double readoutboardSize_X = getDimension(detname+"_X");
  G4double readoutboardSize_Y = getDimension(detname+"_Y");
  G4double readoutboardSize_Z = getDimension(detname+"_Z");
  G4Material* readoutboardMat = G4Material::GetMaterial("PCB");
  G4VisAttributes* readoutboardVisAttr = new G4VisAttributes((nSensorsPerSide == 3 ? G4Colour::Red() : G4Colour::Brown())); readoutboardVisAttr->SetVisibility(true);

  // Connectors/flex/stiffener
  detname = detbase + "_Connectors";
  G4double connectorSize_X = getDimension(detname+"_X");
  G4double connectorSize_Y = getDimension(detname+"_Y");
  G4double connectorSize_Z = getDimension(detname+"_Z");
  G4Material* connectorMat = G4Material::GetMaterial("ServiceConnector");
  G4VisAttributes* connectorVisAttr = new G4VisAttributes((nSensorsPerSide == 3 ? G4Colour::Red() : G4Colour::Brown())); connectorVisAttr->SetVisibility(true);

  // Power board
  detname = detbase + "_PowerBoard";
  G4double powerboardSize_X = getDimension(detname+"_X");
  G4double powerboardSize_Y = getDimension(detname+"_Y");
  G4double powerboardSize_Z = getDimension(detname+"_Z");
  G4Material* powerboardMat = G4Material::GetMaterial("PCB");
  G4VisAttributes* powerboardVisAttr = new G4VisAttributes((nSensorsPerSide == 3 ? G4Colour::Red() : G4Colour::Brown())); powerboardVisAttr->SetVisibility(true);


  // Build the service hybrid
  detname = detbase;
  G4Box* serviceBox = new G4Box(
    detname.c_str(),
    servicehybridSize_X/2., servicehybridSize_Y/2., servicehybridSize_Z/2.
  );
  G4LogicalVolume* serviceLogical = new G4LogicalVolume(
    serviceBox,
    servicehybridMat,
    detname.c_str()
  );
  serviceLogical->SetVisAttributes(servicehybridVisAttr);
  new G4PVPlacement(
    rotation, relativePos,
    serviceLogical,
    detname.c_str(),
    motherLogical,
    false, 0, fCheckOverlaps
  );
  detector_components[detname] = BasicDetectorAttributes(serviceBox, serviceLogical, servicehybridMat, servicehybridVisAttr);


  // Thermal pad
  detname = detbase + "_ThermalPad";
  G4Box* thermalpadBox = new G4Box(
    detname.c_str(),
    thermalpadSize_X/2., thermalpadSize_Y/2., thermalpadSize_Z/2.
  );
  G4LogicalVolume* thermalpadLogical = new G4LogicalVolume(
    thermalpadBox,
    thermalpadMat,
    detname.c_str()
  );
  thermalpadLogical->SetVisAttributes(thermalpadVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, (thermalpadSize_Z-servicehybridSize_Z)/2.),
    thermalpadLogical,
    detname.c_str(),
    serviceLogical,
    false, 0, fCheckOverlaps
  );

  // Readout board
  detname = detbase + "_ReadoutBoard";
  G4Box* readoutboardBox = new G4Box(
    detname.c_str(),
    readoutboardSize_X/2., readoutboardSize_Y/2., readoutboardSize_Z/2.
  );
  G4LogicalVolume* readoutboardLogical = new G4LogicalVolume(
    readoutboardBox,
    readoutboardMat,
    detname.c_str()
  );
  readoutboardLogical->SetVisAttributes(readoutboardVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, thermalpadSize_Z+(readoutboardSize_Z-servicehybridSize_Z)/2.),
    readoutboardLogical,
    detname.c_str(),
    serviceLogical,
    false, 0, fCheckOverlaps
  );

  // Connectors
  detname = detbase + "_Connectors";
  G4Box* connectorBox = new G4Box(
    detname.c_str(),
    connectorSize_X/2., connectorSize_Y/2., connectorSize_Z/2.
  );
  G4LogicalVolume* connectorLogical = new G4LogicalVolume(
    connectorBox,
    connectorMat,
    detname.c_str()
  );
  connectorLogical->SetVisAttributes(connectorVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, thermalpadSize_Z+readoutboardSize_Z+(connectorSize_Z-servicehybridSize_Z)/2.),
    connectorLogical,
    detname.c_str(),
    serviceLogical,
    false, 0, fCheckOverlaps
  );

  // Power board
  detname = detbase + "_PowerBoard";
  G4Box* powerboardBox = new G4Box(
    detname.c_str(),
    powerboardSize_X/2., powerboardSize_Y/2., powerboardSize_Z/2.
  );
  G4LogicalVolume* powerboardLogical = new G4LogicalVolume(
    powerboardBox,
    powerboardMat,
    detname.c_str()
  );
  powerboardLogical->SetVisAttributes(powerboardVisAttr);
  new G4PVPlacement(
    0, G4ThreeVector(0, 0, thermalpadSize_Z+readoutboardSize_Z+connectorSize_Z+(powerboardSize_Z-servicehybridSize_Z)/2.),
    powerboardLogical,
    detname.c_str(),
    serviceLogical,
    false, 0, fCheckOverlaps
  );
}

// The grandmother is the wedge logical, and the mother is the inner passive component
// One needs the grandmother to place the front face support bars.
void DetectorConstruction::BuildWedgeComponents(G4LogicalVolume* passiveLogical, G4LogicalVolume* activeFarLogical, std::vector<std::pair<G4double, G4double>> const& coolingpipes_xpos_ymin){
  using namespace ETLDetectorDimensions;

  string const detbase = "ETLWedge";
  string detname;

  G4Material* empty_mat = G4Material::GetMaterial("Vacuum");
  G4Material* Epoxy_mat = G4Material::GetMaterial("Epoxy");
  G4Material* MIC6Al_mat = G4Material::GetMaterial("MIC6_Al");
  G4Material* CoolingAl_mat = G4Material::GetMaterial("G4_Al");
  G4Material* Cooling_mat = G4Material::GetMaterial("Cool_CO2");
  G4Material* CoolingPipe_mat = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  G4Material* Attachment_mat = CoolingAl_mat;
  if (!Epoxy_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve Epoxy_mat.";
    G4Exception("DetectorConstruction::BuildWedgeComponents",
                "MyCode0001", FatalException, msg);
  }
  if (!MIC6Al_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve MIC6_Al_mat.";
    G4Exception("DetectorConstruction::BuildWedgeComponents",
                "MyCode0001", FatalException, msg);
  }
  if (!CoolingAl_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve CoolingAl_mat.";
    G4Exception("DetectorConstruction::BuildWedgeComponents",
                "MyCode0001", FatalException, msg);
  }
  if (!Cooling_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve Cooling_mat.";
    G4Exception("DetectorConstruction::BuildWedgeComponents",
                "MyCode0001", FatalException, msg);
  }
  if (!CoolingPipe_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve CoolingPipe_mat.";
    G4Exception("DetectorConstruction::BuildWedgeComponents",
                "MyCode0001", FatalException, msg);
  }

  detname = detbase;
  G4double wedge_Rmin = getDimension(detname+"_Rmin");
  G4double wedge_Rmax = getDimension(detname+"_Rmax");
  G4double wedge_MIC6Al_Z = getDimension(detname+"_MIC6Al_Z");
  G4double wedge_Epoxy_Z = getDimension(detname+"_Epoxy_Z");
  G4double wedge_CoolingAl_Z = getDimension(detname+"_CoolingAl_Z");
  G4double wedge_Z = getDimension(detname+"_Z");
  G4double wedge_fullZ = getDimension(detname+"_FullZ");

  detname = detbase + "_Attachment";
  G4double wedge_attachment_X = getDimension(detname+"_X");
  G4double wedge_attachment_Y = getDimension(detname+"_Y");
  G4double wedge_attachment_Z = getDimension(detname+"_Z");
  G4double wedge_attachment_Offset_X = getDimension(detname+"_Offset_X");
  G4double wedge_attachment_Offset_Y = getDimension(detname+"_Offset_Y");
  G4double wedge_attachment_lower_Y = wedge_attachment_Y/2.+wedge_attachment_Offset_Y;
  G4double wedge_attachment_upper_Y = wedge_attachment_Y/2.-wedge_attachment_Offset_Y;

  detname = detbase + "_CoolingPipe";
  G4double coolingpipe_Rmin = getDimension(detname+"_Rmin");
  G4double coolingpipe_Rmax = getDimension(detname+"_Rmax");

  G4double wedge_Rmax_dXOverflow = coolingpipe_Rmax*2.;
  G4double wedge_RSqOverflow = wedge_Rmax_dXOverflow*(wedge_Rmax*2.-wedge_Rmax_dXOverflow);
  G4double wedge_Renclosing = std::sqrt(std::pow(wedge_Rmax, 2)+wedge_RSqOverflow);

  // MIC6 Al
  G4LogicalVolume* logicMIC6AlEnclosing = nullptr;
  if (addMIC6Al){
    detname = detbase + "_MIC6Al";
    G4VisAttributes* MIC6AlBaseVisAttr = new G4VisAttributes(G4Colour::Black()); MIC6AlBaseVisAttr->SetVisibility(false);

    G4Tubs* solidMIC6AlEnclosing = new G4Tubs(
      (detname+"_Enclosing").c_str(),
      wedge_Rmin, wedge_Renclosing, wedge_MIC6Al_Z/2., 0, 90.*degree
    );
    logicMIC6AlEnclosing = new G4LogicalVolume(
      solidMIC6AlEnclosing,
      empty_mat,
      (detname+"_Enclosing").c_str()
    );
    logicMIC6AlEnclosing->SetVisAttributes(MIC6AlBaseVisAttr);
    new G4PVPlacement(
      nullptr,
      G4ThreeVector(0, 0, (-wedge_Z+wedge_MIC6Al_Z)/2.),
      logicMIC6AlEnclosing,
      (detname+"_Enclosing").c_str(),
      passiveLogical,
      false,
      0,
      fCheckOverlaps
    );

    G4Tubs* solidMIC6AlBase = new G4Tubs(
      (detname+"_Base").c_str(),
      wedge_Rmin, wedge_Rmax, wedge_MIC6Al_Z/2., 0, 90.*degree
    );
    G4VSolid* solidMIC6Al = nullptr;
    // The cooling pipes are embedded in the MIC6 Al, so the solid construction involves drilling pipes :(
    if (drillCoolingPipeCavities){
      // Begin construction of the complex solidMIC6Al solid
      std::vector<G4Box*> boxes;
      std::vector<G4Tubs*> tubes;
      std::vector<G4VSolid*> subtractedWedges; subtractedWedges.push_back(solidMIC6AlBase);

      G4double ldy = coolingpipe_Rmax;
      G4double ldx = 2.*ldy;
      G4RotationMatrix rotm; rotm.rotateX(-M_PI/2.*rad);

      for (std::pair<G4double, G4double> const& xpos_ymin:coolingpipes_xpos_ymin){
        G4double const& xpos = xpos_ymin.first;
        G4double const& ymin = xpos_ymin.second;
        //G4double ymax = std::sqrt(std::pow(wedge_Rmax, 2) - std::pow(xpos_ymin.first, 2));
        G4double const& ymax = wedge_Rmax;
        G4double ldz = ymax - ymin;
        G4cout << "Placing a pipe cavity at x = " << xpos_ymin.first << " with y min/max = " << ymin << " / " << ymax << G4endl;
        G4VSolid* subtractionSolid = nullptr;

        //G4Box* theBox = new G4Box((detname+"_ConstituentBox").c_str(), ldx/2., ldy/2., ldz/2.); // dy and dz will be swapped after rotation
        G4Box* theBox = new G4Box((detname+"_ConstituentBox").c_str(), ldx/2., ldx/2., ldz/2.); // dy and dz will be swapped after rotation
        G4cout << "\t- Box dimensions: (x,y,z) = ( " << ldx << ", " << ldx << ", " << ldz << " )" << G4endl;
        boxes.push_back(theBox);
        //G4ThreeVector relpos_box(xpos, ymin+ldz/2., wedge_MIC6Al_Z/2.-ldy/2.);
        G4ThreeVector relpos_box(xpos, ymin+ldz/2., wedge_MIC6Al_Z/2.-ldx/2.);
        G4Transform3D rot_trans_box(rotm, relpos_box);
        subtractionSolid = new G4SubtractionSolid(detname, subtractedWedges.back(), theBox, rot_trans_box);
        subtractedWedges.push_back(subtractionSolid);

        /*
        G4Tubs* theTube = new G4Tubs((detname+"_ConstituentTube").c_str(), G4double(0.), ldy, ldz/2., 0, M_PI*rad*2.);
        G4cout << "\t- Tube dimensions: (r,phi,z) = ( " << ldy << ", " << M_PI*rad*2. << ", " << ldz << " )" << G4endl;
        tubes.push_back(theTube);
        G4ThreeVector relpos_tube(xpos, ymin+ldz/2., wedge_MIC6Al_Z/2.-ldy);
        G4Transform3D rot_trans_tube(rotm, relpos_tube);
        subtractionSolid = new G4SubtractionSolid(detname, subtractedWedges.back(), theTube, rot_trans_tube);
        subtractedWedges.push_back(subtractionSolid);
        */
      }
      solidMIC6Al = subtractedWedges.back();
      // End construction of the complex solidMIC6Al solid
    }
    else{
      // Begin construction of the basic solidMIC6Al solid
      solidMIC6Al = solidMIC6AlBase;
      // End construction of the basic solidMIC6Al solid
    }
    G4LogicalVolume* logicMIC6Al = new G4LogicalVolume(
      solidMIC6Al,
      MIC6Al_mat,
      detname.c_str()
    );
    G4VisAttributes* MIC6AlVisAttr = new G4VisAttributes(G4Colour::Gray()); MIC6AlVisAttr->SetVisibility(true);
    logicMIC6Al->SetVisAttributes(MIC6AlVisAttr);
    new G4PVPlacement(
      nullptr,
      G4ThreeVector(),
      logicMIC6Al,
      detname.c_str(),
      logicMIC6AlEnclosing,
      false,
      0,
      fCheckOverlaps
    );
  }

  if (addCoolingPipes && logicMIC6AlEnclosing){
    detname = detbase + "_CoolingPipe";
    G4RotationMatrix* rotm = new G4RotationMatrix(); rotm->rotateX(M_PI*rad/2.);
    G4VisAttributes* visAttrOuter = new G4VisAttributes(G4Colour::White()); visAttrOuter->SetVisibility(true);
    G4VisAttributes* visAttrInner = new G4VisAttributes(G4Colour::Blue()); visAttrInner->SetVisibility(true);
    for (std::pair<G4double, G4double> const& xpos_ymin:coolingpipes_xpos_ymin){
      G4double const& ymin = xpos_ymin.second;
      G4double ymax = std::sqrt(std::pow(wedge_Rmax, 2) - std::pow(xpos_ymin.first-coolingpipe_Rmax, 2));
      G4double dy = ymax-ymin;
      G4Tubs* solidInnerTube = new G4Tubs((detname+"_Inner").c_str(), G4double(0.), coolingpipe_Rmin, dy/2., G4double(0), M_PI*rad*2.);
      G4Tubs* solidOuterTube = new G4Tubs((detname+"_Outer").c_str(), coolingpipe_Rmin, coolingpipe_Rmax, dy/2., G4double(0), M_PI*rad*2.);
      G4cout << "Placing a pipe at x = " << xpos_ymin.first << " with y min/max = " << ymin << " / " << ymax << " and Rmin/max = " << coolingpipe_Rmin << " / " << coolingpipe_Rmax << G4endl;

      G4LogicalVolume* logicInnerTube = new G4LogicalVolume(
        solidInnerTube,
        Cooling_mat,
        (detname+"_Inner").c_str()
      );
      logicInnerTube->SetVisAttributes(visAttrInner);
      G4LogicalVolume* logicOuterTube = new G4LogicalVolume(
        solidOuterTube,
        CoolingPipe_mat,
        (detname+"_Outer").c_str()
      );
      logicOuterTube->SetVisAttributes(visAttrOuter);

      new G4PVPlacement(
        rotm,
        G4ThreeVector(xpos_ymin.first, ymin+dy/2., wedge_MIC6Al_Z/2.-coolingpipe_Rmax),
        logicInnerTube,
        (detname+"_Inner").c_str(),
        logicMIC6AlEnclosing,
        false,
        0,
        fCheckOverlaps
      );
      new G4PVPlacement(
        rotm,
        G4ThreeVector(xpos_ymin.first, ymin+dy/2., wedge_MIC6Al_Z/2.-coolingpipe_Rmax),
        logicOuterTube,
        (detname+"_Outer").c_str(),
        logicMIC6AlEnclosing,
        false,
        0,
        fCheckOverlaps
      );
    }
  }

  // Epoxy
  detname = detbase + "_Epoxy";
  if (addEpoxy){
    G4Tubs* solidEpoxy = new G4Tubs(
      detname.c_str(),
      wedge_Rmin, wedge_Rmax, wedge_Epoxy_Z/2., 0, 90.*degree
    );
    G4LogicalVolume* logicEpoxy = new G4LogicalVolume(
      solidEpoxy,
      Epoxy_mat,
      detname.c_str()
    );
    G4VisAttributes* EpoxyVisAttr = new G4VisAttributes(G4Colour::Brown()); EpoxyVisAttr->SetVisibility(true);
    logicEpoxy->SetVisAttributes(EpoxyVisAttr);
    new G4PVPlacement(
      nullptr,
      G4ThreeVector(0, 0, (-wedge_Z+wedge_Epoxy_Z)/2.+wedge_MIC6Al_Z),
      logicEpoxy,
      detname.c_str(),
      passiveLogical,
      false,
      0,
      fCheckOverlaps
    );
  }

  // CoolingAl
  detname = detbase + "_CoolingAl";
  if (addCoolingAl){
    G4Tubs* solidCoolingAl = new G4Tubs(
      detname.c_str(),
      wedge_Rmin, wedge_Rmax, wedge_CoolingAl_Z/2., 0, 90.*degree
    );
    G4LogicalVolume* logicCoolingAl = new G4LogicalVolume(
      solidCoolingAl,
      CoolingAl_mat,
      detname.c_str()
    );
    G4VisAttributes* CoolingAlVisAttr = new G4VisAttributes(G4Colour::Gray()); CoolingAlVisAttr->SetVisibility(true);
    logicCoolingAl->SetVisAttributes(CoolingAlVisAttr);
    new G4PVPlacement(
      nullptr,
      G4ThreeVector(0, 0, (-wedge_Z+wedge_CoolingAl_Z)/2.+wedge_MIC6Al_Z+wedge_Epoxy_Z),
      logicCoolingAl,
      detname.c_str(),
      passiveLogical,
      false,
      0,
      fCheckOverlaps
    );
  }

  detname = detbase + "_Attachment";
  if (addFrontFaceSupportBars){
    G4VisAttributes* AttachmentVisAttr = new G4VisAttributes(G4Colour::Gray()); AttachmentVisAttr->SetVisibility(true);

    G4ThreeVector relpos_lower(wedge_attachment_Offset_X + wedge_attachment_X/2., wedge_attachment_lower_Y/2., -(wedge_fullZ-wedge_Z)/4.+wedge_attachment_Z/2.);
    G4Box* solidAttachment_Lower = new G4Box((detname+"_Lower").c_str(), wedge_attachment_X/2., wedge_attachment_lower_Y/2., wedge_attachment_Z/2.);
    G4LogicalVolume* logicAttachment_Lower = new G4LogicalVolume(
      solidAttachment_Lower,
      Attachment_mat,
      (detname+"_Lower").c_str()
    );
    logicAttachment_Lower->SetVisAttributes(AttachmentVisAttr);
    new G4PVPlacement(
      nullptr,
      relpos_lower,
      logicAttachment_Lower,
      (detname+"_Lower").c_str(),
      activeFarLogical,
      false,
      0,
      fCheckOverlaps
    );

    G4ThreeVector relpos_upper(wedge_attachment_upper_Y/2., wedge_attachment_Offset_X + wedge_attachment_X/2., -(wedge_fullZ-wedge_Z)/4.+wedge_attachment_Z/2.);
    G4RotationMatrix* rotm = new G4RotationMatrix; rotm->rotateZ(-M_PI*rad/2.);
    G4Box* solidAttachment_Upper = new G4Box((detname+"_Upper").c_str(), wedge_attachment_X/2., wedge_attachment_upper_Y/2., wedge_attachment_Z/2.);
    G4LogicalVolume* logicAttachment_Upper = new G4LogicalVolume(
      solidAttachment_Upper,
      Attachment_mat,
      (detname+"_Upper").c_str()
    );
    logicAttachment_Upper->SetVisAttributes(AttachmentVisAttr);
    new G4PVPlacement(
      rotm,
      relpos_upper,
      logicAttachment_Upper,
      (detname+"_Upper").c_str(),
      activeFarLogical,
      false,
      0,
      fCheckOverlaps
    );
  }
}

void DetectorConstruction::BuildEndCapSupportComponents(G4LogicalVolume* motherLogical){
  using namespace ETLDetectorDimensions;

  string const detbase = "ETLSupport";
  string const det_wedge = "ETLWedge";
  string const det_offset = "ETLOffset";
  string detname;

  G4Material* empty_mat = G4Material::GetMaterial("Vacuum");
  G4Material* supportTube_mat = G4Material::GetMaterial("Permaglas");
  G4Material* neutronModerator_mat = G4Material::GetMaterial("G4_POLYETHYLENE");
  G4Material* thermalScreenCore_mat = G4Material::GetMaterial("Aerogel");
  G4Material* thermalScreenSkin_mat = supportTube_mat;
  G4Material* mountingBracket_mat = supportTube_mat;
  if (!empty_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve empty_mat.";
    G4Exception("DetectorConstruction::BuildEndCapSupportComponents",
                "MyCode0001", FatalException, msg);
  }
  if (!supportTube_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve supportTube_mat.";
    G4Exception("DetectorConstruction::BuildEndCapSupportComponents",
                "MyCode0001", FatalException, msg);
  }
  if (!neutronModerator_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve neutronModerator_mat.";
    G4Exception("DetectorConstruction::BuildEndCapSupportComponents",
                "MyCode0001", FatalException, msg);
  }
  if (!thermalScreenCore_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve thermalScreenCore_mat.";
    G4Exception("DetectorConstruction::BuildEndCapSupportComponents",
                "MyCode0001", FatalException, msg);
  }

  detname = det_wedge;
  G4double wedge_Rmin = getDimension(detname+"_Rmin");
  G4double wedge_Z = getDimension(detname+"_Z");
  G4double wedge_fullZ = getDimension(detname+"_FullZ");

  detname = det_offset + "_Disks_dZ";
  G4double sep_Z_disks = getDimension(detname);

  detname = detbase + "_SupportTube";
  G4double supportTube_thickness = getDimension(detname+"_Thickness");
  G4double const& supportTube_Rmax = wedge_Rmin;
  G4double supportTube_Rmin = supportTube_Rmax - supportTube_thickness;
  G4double supportTube_Z = sep_Z_disks + wedge_Z + wedge_fullZ;

  detname = detbase + "_MountingBracket";
  G4double mountingBracket_thickness = getDimension(detname+"_Thickness");
  G4double const& mountingBracket_Rmin = wedge_Rmin;
  G4double mountingBracket_Rmax = mountingBracket_Rmin + mountingBracket_thickness;

  // Support tube
  detname = detbase + "_SupportTube";
  if (addSupportTube){
    G4VisAttributes* supportTubeVisAttr = new G4VisAttributes(G4Colour::Blue()); supportTubeVisAttr->SetVisibility(true);
    G4Tubs* solidSupportTube = new G4Tubs(
      detname.c_str(),
      supportTube_Rmin, supportTube_Rmax, supportTube_Z/2., 0, M_PI*2.*rad
    );
    G4LogicalVolume* logicSupportTube = new G4LogicalVolume(
      solidSupportTube,
      supportTube_mat,
      detname.c_str()
    );
    logicSupportTube->SetVisAttributes(supportTubeVisAttr);
    new G4PVPlacement(
      nullptr,
      G4ThreeVector(),
      logicSupportTube,
      detname.c_str(),
      motherLogical,
      false,
      0,
      fCheckOverlaps
    );
  }

  // Mounting bracket on top of the support tube and in the middle of the two disks
  detname = detbase + "_MountingBracket";
  if (addMountingBracket){
    G4VisAttributes* mountingBracketVisAttr = new G4VisAttributes(G4Colour::Blue()); mountingBracketVisAttr->SetVisibility(true);
    G4Tubs* solidMountingBracket = new G4Tubs(
      detname.c_str(),
      mountingBracket_Rmin, mountingBracket_Rmax, sep_Z_disks/2., 0, M_PI*2.*rad
    );
    G4LogicalVolume* logicMountingBracket = new G4LogicalVolume(
      solidMountingBracket,
      mountingBracket_mat,
      detname.c_str()
    );
    logicMountingBracket->SetVisAttributes(mountingBracketVisAttr);
    new G4PVPlacement(
      nullptr,
      G4ThreeVector(),
      logicMountingBracket,
      detname.c_str(),
      motherLogical,
      false,
      0,
      fCheckOverlaps
    );
  }

}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes(){
  using namespace ETLDetectorDimensions;

  // Setup materials
  G4Material* world_mat = G4Material::GetMaterial("Vacuum");
  G4Material* diskBox_mat = world_mat;
  G4Material* wedgeEnclosure_mat = diskBox_mat;
  G4Material* endcap_mat = diskBox_mat;
  G4Material* mtd_mat = diskBox_mat;
  G4Material* wedge_attachment_mat = G4Material::GetMaterial("G4_Al");
  if (!world_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve world_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }
  if (!diskBox_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve diskBox_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }
  if (!wedge_attachment_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve wedge_attachment_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }

  // Size parameters
  string detname;
  string const det_cms = "CMS";
  string const det_mtd = "MTD";
  string const det_etl = "MTD_ETL";
  string const det_endcap = "ETLEndCap";
  string const det_disk = "ETLDisk";
  string const det_halfdisk = "ETLHalfDisk";
  string const det_offset = "ETLOffset";
  string const det_wedge = "ETLWedge";
  string const det_wedge_attachment = "ETLWedge_Attachment";
  string const det_wedge_coolingpipe = "ETLWedge_CoolingPipe";
  string const det_onesensor = "ETLOneSensorModule";
  string const det_twosensor = "ETLTwoSensorModule";
  string const det_servicehybrid6 = "ETL6SensorServiceHybrid";
  string const det_servicehybrid12 = "ETL12SensorServiceHybrid";
  string const det_support = "ETLSupport";

  detname = det_onesensor;
  G4double onesensor_X = getDimension(detname+"_X");
  G4double onesensor_Y = getDimension(detname+"_Y");
  G4double onesensor_Z = getDimension(detname+"_Z");

  detname = det_twosensor;
  G4double twosensor_X = getDimension(detname+"_X");
  G4double twosensor_Y = getDimension(detname+"_Y");
  G4double twosensor_Z = getDimension(detname+"_Z");
  if (onesensor_Y != twosensor_Y){
    G4ExceptionDescription msg;
    msg << "Inconsistent one- and two-sensor y-dimensions!";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0002", FatalException, msg);
  }

  detname = det_servicehybrid6;
  G4double servicehybrid6_X = getDimension(detname+"_X");
  G4double servicehybrid6_Y = getDimension(detname+"_Y");
  G4double servicehybrid6_Z = getDimension(detname+"_Z");

  detname = det_servicehybrid12;
  G4double servicehybrid12_X = getDimension(detname+"_X");
  G4double servicehybrid12_Y = getDimension(detname+"_Y");
  G4double servicehybrid12_Z = getDimension(detname+"_Z");
  if (servicehybrid6_X != servicehybrid12_X){
    G4ExceptionDescription msg;
    msg << "Inconsistent service hybrid x-dimensions!";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0002", FatalException, msg);
  }

  detname = det_wedge;
  G4double wedge_Rmin = getDimension(detname+"_Rmin");
  G4double wedge_Rmax = getDimension(detname+"_Rmax");
  G4double wedge_Z = getDimension(detname+"_Z");
  G4double wedge_fullZ = getDimension(detname+"_FullZ");

  // Thin aluminum attachment on the wedge front faces
  detname = det_wedge_attachment;
  G4double wedge_attachment_X = getDimension(detname+"_X");
  G4double wedge_attachment_Y = getDimension(detname+"_Y");
  G4double wedge_attachment_Z = getDimension(detname+"_Z");
  G4double wedge_attachment_Offset_X = getDimension(detname+"_Offset_X");
  G4double wedge_attachment_Offset_Y = getDimension(detname+"_Offset_Y");
  G4double wedge_frontface_sensor_Offset_X = wedge_attachment_Y/2. - wedge_attachment_Offset_Y;
  G4double wedge_frontface_sensor_Offset_Y = wedge_attachment_Y/2. + wedge_attachment_Offset_Y;

  // Cooling pipes
  detname = det_wedge_coolingpipe;
  G4double coolingpipe_Rmax = getDimension(detname+"_Rmax");

  // Support structures
  detname = det_support;
  G4double mountingBracket_thickness = getDimension(detname+"_MountingBracket_Thickness");

  // External offsets
  detname = det_offset + "_Module_SensorServiceHybrid_dX";
  G4double sep_X_module_servicehybrid = getDimension(detname);
  detname = det_offset + "_Module_Module_dY";
  G4double sep_Y_module_module = getDimension(detname);
  detname = det_offset + "_Module_WedgeAttachment_dX";
  G4double sep_X_module_wedgeAttachment = getDimension(detname);
  detname = det_offset + "_Module_WedgeAttachment_dY";
  G4double sep_Y_module_wedgeAttachment = getDimension(detname);
  detname = det_offset + "_SensorServiceHybrid_SensorServiceHybrid_dY";
  G4double sep_Y_servicehybrid_servicehybrid = getDimension(detname);
  detname = det_offset + "_Disks_dZ";
  G4double sep_Z_disks = getDimension(detname);
  detname = det_offset + "_IP_dZ";
  G4double IP_dZ = getDimension(detname);


  /******************/
  /******************/
  /* BEGIN GEOMETRY */
  /******************/
  /******************/

  // Generic wedge box parameters
  G4double wedge_Rmax_dXOverflow = coolingpipe_Rmax*2.;
  G4double wedge_RSqOverflow = std::max(wedge_Rmax_dXOverflow*(wedge_Rmax*2.-wedge_Rmax_dXOverflow), std::pow(wedge_attachment_Y, 2));
  G4double wedge_Renclosing = std::sqrt(std::pow(wedge_Rmax, 2)+wedge_RSqOverflow);
  G4double const& diskBox_Rmin = wedge_Rmin;
  G4double const& diskBox_Rmax = wedge_Renclosing;
  G4double diskBox_X = wedge_Renclosing*2.;
  G4double diskBox_Y = wedge_Renclosing*2.;
  G4double diskBox_Z = wedge_fullZ;

  // Generic disk pair parameters
  G4double diskPair_X = diskBox_X;
  G4double diskPair_Y = diskBox_Y;
  G4double diskPair_Z = diskBox_Z*2. + sep_Z_disks;

  // Generic endcap parameters
  G4double endcap_X = diskPair_X;
  G4double endcap_Y = diskPair_Y;
  G4double endcap_Z = diskPair_Z; // Will have to change once other systems are also placed

  // Generic MTD parameters
  G4double mtd_X = endcap_X;
  G4double mtd_Y = endcap_Y;
  G4double mtd_Z = diskPair_Z*2. + IP_dZ + (endcap_Z-diskPair_Z)*2.;

  // World parameters
  G4double world_X = mtd_X;
  G4double world_Y = mtd_Y;
  G4double world_Z = mtd_Z;

  // World
  detname = det_cms;
  G4Box* solidWorld = new G4Box(
    detname.c_str(),
    world_X/2., world_Y/2., world_Z/2.
  );
  G4LogicalVolume* logicWorld = new G4LogicalVolume(
    solidWorld,
    world_mat,
    detname.c_str()
  );
  G4VisAttributes* worldVisAttr = new G4VisAttributes(G4Colour::Black()); worldVisAttr->SetVisibility(false);
  logicWorld->SetVisAttributes(worldVisAttr);
  G4VPhysicalVolume* physWorld = new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logicWorld,
    detname.c_str(),
    nullptr,
    false,
    0,
    fCheckOverlaps
  );

  // The box enclosure of the two ETLs
  detname = det_mtd;
  G4Box* solidMTD = new G4Box(
    detname.c_str(),
    mtd_X/2., mtd_Y/2., mtd_Z/2.
  );
  G4LogicalVolume* logicMTD = new G4LogicalVolume(
    solidMTD,
    mtd_mat,
    detname.c_str()
  );
  G4VisAttributes* mtdVisAttr = new G4VisAttributes(G4Colour::Black()); mtdVisAttr->SetVisibility(false);
  logicMTD->SetVisAttributes(mtdVisAttr);
  new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logicMTD,
    detname.c_str(),
    logicWorld,
    false,
    0,
    fCheckOverlaps
  );
  detector_components[detname] = BasicDetectorAttributes(solidMTD, logicMTD, mtd_mat, mtdVisAttr);

  // The box enclosure of ETL
  detname = det_endcap;
  G4Box* solidEndcap = new G4Box(
    detname.c_str(),
    endcap_X/2., endcap_Y/2., endcap_Z/2.
  );
  G4LogicalVolume* logicEndcap = new G4LogicalVolume(
    solidEndcap,
    endcap_mat,
    detname.c_str()
  );
  G4VisAttributes* endcapVisAttr = new G4VisAttributes(G4Colour::Black()); endcapVisAttr->SetVisibility(false);
  logicEndcap->SetVisAttributes(endcapVisAttr);
  detector_components[detname] = BasicDetectorAttributes(solidEndcap, logicEndcap, endcap_mat, endcapVisAttr);
  // Defer placement until the very end

  // The box enclosure of a full disk
  detname = det_disk;
  G4Tubs* solidDiskBox = new G4Tubs(
    detname.c_str(),
    diskBox_Rmin, diskBox_Rmax, diskBox_Z/2., 0, M_PI*2.*rad
  );
  G4LogicalVolume* logicDiskBox = new G4LogicalVolume(
    solidDiskBox,
    diskBox_mat,
    detname.c_str()
  );
  G4VisAttributes* diskBoxVisAttr = new G4VisAttributes(G4Colour::Black()); diskBoxVisAttr->SetVisibility(false);
  logicDiskBox->SetVisAttributes(diskBoxVisAttr);
  detector_components[detname] = BasicDetectorAttributes(solidDiskBox, logicDiskBox, diskBox_mat, diskBoxVisAttr);
  // Defer the placement of the disk boxes until the end

  detname = det_disk+"_WedgeContainer";
  G4Tubs* solidDisk = new G4Tubs(
    detname.c_str(),
    wedge_Rmin, wedge_Renclosing, wedge_fullZ/2., 0, M_PI*2.*rad
  );
  G4LogicalVolume* logicDisk = new G4LogicalVolume(
    solidDisk,
    diskBox_mat,
    detname.c_str()
  );
  logicDisk->SetVisAttributes(diskBoxVisAttr);
  detector_components[detname] = BasicDetectorAttributes(solidDisk, logicDisk, diskBox_mat, diskBoxVisAttr);
  new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logicDisk,
    detname.c_str(),
    logicDiskBox,
    false,
    0,
    fCheckOverlaps
  );

  // Wedge
  detname = det_wedge;
  G4Tubs* solidWedge = new G4Tubs(
    detname.c_str(),
    wedge_Rmin, wedge_Renclosing, wedge_fullZ/2., 0, 90.*degree
  );
  G4LogicalVolume* logicWedge = new G4LogicalVolume(
    solidWedge,
    wedgeEnclosure_mat, // Notice here that the material defined is the dummy disk box material, not the actual materials for the wedge components
    detname.c_str()
  );
  G4VisAttributes* wedgeVisAttr = new G4VisAttributes(G4Colour::Black()); wedgeVisAttr->SetVisibility(false); // The wedge object should not be visible
  logicWedge->SetVisAttributes(wedgeVisAttr);
  if (!doFullDisk) new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logicWedge,
    detname.c_str(),
    logicDisk,
    false,
    0,
    fCheckOverlaps
  );
  else new G4PVReplica(
    detname.c_str(),
    logicWedge,
    logicDisk,
    kPhi, 4, M_PI*rad/2., -M_PI*rad/4.
  );
  detector_components[detname] = BasicDetectorAttributes(solidWedge, logicWedge, wedgeEnclosure_mat, wedgeVisAttr);

  // Wedge passive material
  detname = det_wedge+"_PassiveMaterial";
  G4Tubs* solidWedgePassive = new G4Tubs(
    detname.c_str(),
    wedge_Rmin, wedge_Renclosing, wedge_Z/2., 0, 90.*degree
  );
  G4LogicalVolume* logicWedgePassive = new G4LogicalVolume(
    solidWedgePassive,
    wedgeEnclosure_mat, // Notice here that the material defined is the dummy disk box material, not the actual materials for the wedge components
    detname.c_str()
  );
  logicWedgePassive->SetVisAttributes(wedgeVisAttr);
  detector_components[detname] = BasicDetectorAttributes(solidWedgePassive, logicWedgePassive, wedgeEnclosure_mat, wedgeVisAttr);
  // Defer the placement of the passive material

  // Wedge active material
  detname = det_wedge+"_ActiveMaterial";
  G4Tubs* solidWedgeActive = new G4Tubs(
    detname.c_str(),
    wedge_Rmin, wedge_Renclosing, (wedge_fullZ-wedge_Z)/4., 0, 90.*degree
  );
  G4LogicalVolume* logicWedgeActive_Far = new G4LogicalVolume(
    solidWedgeActive,
    wedgeEnclosure_mat, // Notice here that the material defined is the dummy disk box material, not the actual materials for the wedge components
    (detname+"_Far").c_str()
  );
  logicWedgeActive_Far->SetVisAttributes(wedgeVisAttr);
  detector_components[(detname+"_Far")] = BasicDetectorAttributes(solidWedgeActive, logicWedgeActive_Far, wedgeEnclosure_mat, wedgeVisAttr);
  G4LogicalVolume* logicWedgeActive_Close = new G4LogicalVolume(
    solidWedgeActive,
    wedgeEnclosure_mat, // Notice here that the material defined is the dummy disk box material, not the actual materials for the wedge components
    (detname+"_Close").c_str()
  );
  logicWedgeActive_Close->SetVisAttributes(wedgeVisAttr);
  detector_components[(detname+"_Close")] = BasicDetectorAttributes(solidWedgeActive, logicWedgeActive_Close, wedgeEnclosure_mat, wedgeVisAttr);
  // Defer the placement of the active material

  // Variables for the construction of service hybrids and the modules
  size_t ix_module, ix_service;
  G4double wedge_xpos;
  G4double wedge_ypos;
  G4double wedge_yposinf;
  G4double wedge_yposmin;
  G4double wedge_yposmax;
  G4double wedge_Roffset;
  G4RotationMatrix* reflectionTransformation;
  G4RotationMatrix* reflectionAndSideSwapTransformation = new G4RotationMatrix;
  std::vector<std::pair<G4double, G4double>> coolingpipes_xpos_ymin;
  std::vector<std::pair<G4double, G4double>> sensorhybrid_yminmax;
  std::vector<G4double> sensorhybrid_xpos;
  std::vector<std::pair<G4double, G4double>>::const_iterator moduleSensorHybridConnection_yminmax_left, moduleSensorHybridConnection_yminmax_right, moduleSensorHybridConnection_yminmax_cend;
  std::vector<G4double>::const_iterator moduleSensorHybridConnection_xpos_left, moduleSensorHybridConnection_xpos_right, moduleSensorHybridConnection_xpos_cend;


  /*****************************************/
  /*****************************************/
  /* Construct the front face of the wedge */
  /*****************************************/
  /*****************************************/
  wedge_Roffset = 0;
  ix_module=ix_service=0;
  reflectionTransformation = nullptr;
  reflectionAndSideSwapTransformation = new G4RotationMatrix; reflectionAndSideSwapTransformation->rotateZ(M_PI*rad);

  // Place service hybrids first
  wedge_xpos = wedge_frontface_sensor_Offset_X + sep_X_module_wedgeAttachment + onesensor_X + sep_X_module_servicehybrid + servicehybrid6_X/2.;
  // Calculate yposmin:
  // Depends on how many modules will fit as well
  wedge_yposmin = (wedge_Rmin + wedge_Roffset);
  wedge_yposinf = (wedge_frontface_sensor_Offset_Y + sep_Y_module_wedgeAttachment);
  {
    G4double nmodules_d = (wedge_yposmin - wedge_yposinf) / (onesensor_Y + sep_Y_module_module);
    G4int nmodules_i = std::ceil(nmodules_d);
    wedge_yposmin = ((G4double) nmodules_i) * (onesensor_Y + sep_Y_module_module) + wedge_yposinf;
  }
  while (doWedgeFarFace){ // Loop over columns
    sensorhybrid_xpos.push_back(wedge_xpos);
    coolingpipes_xpos_ymin.emplace_back(wedge_xpos, wedge_yposmin);

    wedge_yposmax = sqrt(fabs(pow(wedge_Rmax, 2) - pow(wedge_xpos + servicehybrid6_X/2., 2)));
    G4cout << "y min/max = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;
    sensorhybrid_yminmax.push_back(std::pair<G4double, G4double>(wedge_yposmin, wedge_yposmin)); // wedge_yposmax is not exactly the location of the edge.
    std::pair<G4double, G4double>& last_sensorhybrid_yminmax = sensorhybrid_yminmax.back();

    size_t i_object=0;
    wedge_ypos = wedge_yposmin;
    while (true){
      size_t n_modules_per_side;
      if (
        (ix_service<5 && i_object==0)
        ||
        (ix_service==3 && i_object<3)
        ||
        ((wedge_yposmax-wedge_ypos)<servicehybrid12_Y && (wedge_yposmax-wedge_ypos)>=servicehybrid6_Y)
        ) n_modules_per_side = 3;
      else n_modules_per_side = 6;

      bool doPlaceSensorHybrid=true;

      switch (n_modules_per_side){
      case 3:
      {
        wedge_ypos += servicehybrid6_Y/2.;
        if (wedge_ypos + servicehybrid6_Y/2.>=wedge_yposmax){ doPlaceSensorHybrid = false; break; } // Breaks from the switch, not the while loop!

        // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, -(wedge_fullZ-wedge_Z)/4. + servicehybrid6_Z/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid6_X << "," << servicehybrid6_Y << "," << servicehybrid6_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 6-sensor service hybrid
        if (putServiceHybrids) BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedgeActive_Far, reflectionTransformation, relpos
        );

        wedge_ypos += servicehybrid6_Y/2.;
        last_sensorhybrid_yminmax.second += servicehybrid6_Y;
        break;
      }
      case 6:
      {
        wedge_ypos += servicehybrid12_Y/2.;
        if (wedge_ypos + servicehybrid12_Y/2.>=wedge_yposmax){ doPlaceSensorHybrid = false; break; } // Breaks from the switch, not the while loop!

        // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, -(wedge_fullZ-wedge_Z)/4. + servicehybrid12_Z/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid12_X << "," << servicehybrid12_Y << "," << servicehybrid12_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 12-sensor service hybrid
        if (putServiceHybrids) BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedgeActive_Far, reflectionTransformation, relpos
        );

        wedge_ypos += servicehybrid12_Y/2.;
        last_sensorhybrid_yminmax.second += servicehybrid12_Y;
        break;
      }
      default:
      {
        G4ExceptionDescription msg;
        msg << "n_modules_per_side = " << n_modules_per_side << " is not supported!";
        G4Exception("DetectorConstruction::DefineVolumes",
                    "MyCode0003", FatalException, msg);
        break;
      }
      }
      i_object++;

      wedge_ypos += sep_Y_servicehybrid_servicehybrid;
      if (doPlaceSensorHybrid) last_sensorhybrid_yminmax.second += sep_Y_servicehybrid_servicehybrid;
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_service++;
    last_sensorhybrid_yminmax.second -= sep_Y_servicehybrid_servicehybrid;

    wedge_xpos += (twosensor_X + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    while (
      wedge_yposmin>0.
      && (
        (wedge_xpos - servicehybrid6_X/2.)>=(wedge_Rmin + wedge_Roffset)
        ||
        (wedge_yposmin - sqrt(fabs(pow(wedge_Rmin + wedge_Roffset, 2) - pow(wedge_xpos - servicehybrid6_X/2., 2))))>=(onesensor_Y+sep_Y_module_module)
        )
      ){
      if (wedge_yposmin<(onesensor_Y + sep_Y_module_module)) break;
      G4cout << "Subtracting deltaY = " << (onesensor_Y+sep_Y_module_module) << " from y pos = " << wedge_yposmin << G4endl;
      wedge_yposmin -= (onesensor_Y + sep_Y_module_module);
      G4cout << "New y pos = " << wedge_yposmin << G4endl;
    }

    if (wedge_xpos>=wedge_Rmax) break;
  }

  // Place one or two-sensor modules next
  wedge_xpos = wedge_frontface_sensor_Offset_X + sep_X_module_wedgeAttachment + onesensor_X/2.;
  moduleSensorHybridConnection_yminmax_cend = sensorhybrid_yminmax.cend();
  moduleSensorHybridConnection_yminmax_left = moduleSensorHybridConnection_yminmax_cend;
  moduleSensorHybridConnection_yminmax_right = sensorhybrid_yminmax.cbegin();
  moduleSensorHybridConnection_xpos_cend = sensorhybrid_xpos.cend();
  moduleSensorHybridConnection_xpos_left = moduleSensorHybridConnection_xpos_cend;
  moduleSensorHybridConnection_xpos_right = sensorhybrid_xpos.cbegin();
  while (doWedgeFarFace){ // Loop over columns
    if (moduleSensorHybridConnection_yminmax_left!=moduleSensorHybridConnection_yminmax_cend && moduleSensorHybridConnection_yminmax_right!=moduleSensorHybridConnection_yminmax_cend){
      wedge_yposmin = std::min(moduleSensorHybridConnection_yminmax_left->first, moduleSensorHybridConnection_yminmax_right->first);
      wedge_yposmax = std::max(moduleSensorHybridConnection_yminmax_left->second, moduleSensorHybridConnection_yminmax_right->second);
    }
    else if (moduleSensorHybridConnection_yminmax_left!=moduleSensorHybridConnection_yminmax_cend){
      wedge_yposmin = moduleSensorHybridConnection_yminmax_left->first;
      wedge_yposmax = moduleSensorHybridConnection_yminmax_left->second;
    }
    else if (moduleSensorHybridConnection_yminmax_right!=moduleSensorHybridConnection_yminmax_cend){
      wedge_yposmin = moduleSensorHybridConnection_yminmax_right->first;
      wedge_yposmax = moduleSensorHybridConnection_yminmax_right->second;
    }
    G4cout << "Sensor y inf/sup = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;

    bool firstColumn = (wedge_xpos<sensorhybrid_xpos.front());
    bool firstColumn_OneSensorModule = (firstColumn && sensorhybrid_xpos.front()>=onesensor_X+sep_X_module_servicehybrid);
    bool lastColumn = (wedge_xpos>sensorhybrid_xpos.back());
    bool lastColumn_OneSensorModule = (lastColumn && (wedge_Rmax-sensorhybrid_xpos.back())>=onesensor_X+sep_X_module_servicehybrid);
    bool useOneSensorModule = (firstColumn_OneSensorModule || lastColumn_OneSensorModule);
    G4double const& moduleWidthX = (useOneSensorModule ? onesensor_X : twosensor_X);
    G4double const& moduleWidthY = (useOneSensorModule ? onesensor_Y : twosensor_Y);
    G4cout << "Placing a set of " << (useOneSensorModule ? "one-" : "two-") << "sensor modules with width " << moduleWidthX << G4endl;

    coolingpipes_xpos_ymin.emplace_back(wedge_xpos, wedge_yposmin);
    // Correct ymin for the common segment of service hybrids
    if (moduleSensorHybridConnection_yminmax_left!=moduleSensorHybridConnection_yminmax_cend && moduleSensorHybridConnection_yminmax_right!=moduleSensorHybridConnection_yminmax_cend)
      coolingpipes_xpos_ymin.back().second = std::max(moduleSensorHybridConnection_yminmax_left->first, moduleSensorHybridConnection_yminmax_right->first);

    size_t i_object=0;
    wedge_ypos = wedge_yposmin;
    while (true){
      G4double wedge_xpos_offset = 0; // This is an offset to calculate if placing left- or right-flanked one-sensor modules in an otherwise two-sensor module column
      bool placeLeftFlankedOneSensorModule = false;
      bool placeRightFlankedOneSensorModule = false;
      bool placeRegularModule = true;
      if (!useOneSensorModule){
        // Check if a left- or right-flanked one-sensor modules needs to be placed. If so, turn off the regular module flag.
        if (
          moduleSensorHybridConnection_yminmax_left != moduleSensorHybridConnection_yminmax_cend
          &&
          moduleSensorHybridConnection_yminmax_right != moduleSensorHybridConnection_yminmax_cend
          ){
          placeRightFlankedOneSensorModule = (
            (wedge_ypos>=moduleSensorHybridConnection_yminmax_left->second && wedge_ypos<moduleSensorHybridConnection_yminmax_right->second)
            ||
            (wedge_ypos<moduleSensorHybridConnection_yminmax_left->first && wedge_ypos>=moduleSensorHybridConnection_yminmax_right->first)
            );
          placeLeftFlankedOneSensorModule = (
            (wedge_ypos<moduleSensorHybridConnection_yminmax_left->second && wedge_ypos>=moduleSensorHybridConnection_yminmax_right->second)
            ||
            (wedge_ypos>=moduleSensorHybridConnection_yminmax_left->first && wedge_ypos<moduleSensorHybridConnection_yminmax_right->first)
            );
          placeRegularModule = !(placeRightFlankedOneSensorModule || placeLeftFlankedOneSensorModule);
        }
        // Check if the placement of a left- or right-flanekd one-sensor module makes sense.
        // Notice placeRegularModule is off if these are true.
        if (
          placeRightFlankedOneSensorModule
          &&
          std::pow(wedge_ypos, 2) + std::pow(wedge_xpos + moduleWidthX/2. - onesensor_X, 2) < pow(wedge_Rmin + wedge_Roffset, 2)
          ) placeRightFlankedOneSensorModule = false;
        if (
          placeLeftFlankedOneSensorModule
          &&
          std::pow(wedge_ypos + onesensor_Y, 2) + std::pow(wedge_xpos - moduleWidthX/2. + onesensor_X, 2) > pow(wedge_Rmax, 2)
          ) placeLeftFlankedOneSensorModule = false;
      }

      if (
        placeRegularModule
        && (
          std::pow(wedge_ypos, 2) + std::pow(wedge_xpos - moduleWidthX/2., 2) < pow(wedge_Rmin + wedge_Roffset, 2)
          ||
          std::pow(wedge_ypos + moduleWidthY, 2) + std::pow(wedge_xpos + moduleWidthX/2., 2) > pow(wedge_Rmax, 2)
          )
        ) placeRegularModule = false;

      if (placeLeftFlankedOneSensorModule) wedge_xpos_offset = (-moduleWidthX+onesensor_X)/2.;
      else if (placeRightFlankedOneSensorModule) wedge_xpos_offset = (+moduleWidthX-onesensor_X)/2.;
      if (placeLeftFlankedOneSensorModule)
        G4cout << "\t- Special left-flanked placement with an x-offset of " << wedge_xpos_offset << " relative to the column center." << G4endl;
      else if (placeRightFlankedOneSensorModule)
        G4cout << "\t- Special right-flanked placement with an x-offset of " << wedge_xpos_offset << " relative to the column center." << G4endl;
      else if (placeRegularModule)
        G4cout << "\t- No special placement; x-offset is " << wedge_xpos_offset << " relative to the column center." << G4endl;

      if ((useOneSensorModule && placeRegularModule) || placeLeftFlankedOneSensorModule || placeRightFlankedOneSensorModule){ // One-sensor modules
        wedge_ypos += onesensor_Y/2.;
        if (wedge_ypos + onesensor_Y/2.>wedge_yposmax) break; // Breaks from the switch, not the while loop!

        // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos+wedge_xpos_offset, wedge_ypos, -(wedge_fullZ-wedge_Z)/4. + onesensor_Z/2.);
        G4cout << "\t- Placing one sensor modules of sizes (" << onesensor_X << "," << onesensor_Y << "," << onesensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        if (putModules) BuildOneSensorModule(
          logicWedgeActive_Far, ((firstColumn_OneSensorModule || placeRightFlankedOneSensorModule) ? reflectionTransformation : reflectionAndSideSwapTransformation), relpos
        );

        wedge_ypos += onesensor_Y/2.;
        i_object++;
      }
      else if (!useOneSensorModule && placeRegularModule){ // Two-sensor modules
        wedge_ypos += twosensor_Y/2.;
        if (wedge_ypos + twosensor_Y/2.>wedge_yposmax) break; // Breaks from the switch, not the while loop!

        // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos+wedge_xpos_offset, wedge_ypos, -(wedge_fullZ-wedge_Z)/4. + twosensor_Z/2.);
        G4cout << "\t- Placing two sensor modules of sizes (" << twosensor_X << "," << twosensor_Y << "," << twosensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        if (putModules) BuildTwoSensorModule(
          logicWedgeActive_Far, reflectionTransformation, relpos
        );

        wedge_ypos += twosensor_Y/2.;
        i_object++;
      }
      else wedge_ypos += (useOneSensorModule ? onesensor_Y : twosensor_Y);

      wedge_ypos += sep_Y_module_module;
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_module++;

    wedge_xpos += (moduleWidthX/2. + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    G4double const& next_moduleWidthX = (wedge_xpos<sensorhybrid_xpos.back() ? twosensor_X : onesensor_X);
    wedge_xpos += next_moduleWidthX/2.;

    if (wedge_xpos+next_moduleWidthX/2.>=wedge_Rmax) break;
    if (moduleSensorHybridConnection_yminmax_right==moduleSensorHybridConnection_yminmax_cend) break;

    moduleSensorHybridConnection_yminmax_left = moduleSensorHybridConnection_yminmax_right;
    moduleSensorHybridConnection_yminmax_right++;
    moduleSensorHybridConnection_xpos_left = moduleSensorHybridConnection_xpos_right;
    moduleSensorHybridConnection_xpos_right++;
  }


  /****************************************/
  /****************************************/
  /* Construct the back face of the wedge */
  /****************************************/
  /****************************************/
  wedge_Roffset = mountingBracket_thickness;
  ix_module=ix_service=0;
  wedge_xpos=0;
  reflectionTransformation = new G4RotationMatrix; reflectionTransformation->rotateY(M_PI*rad);
  reflectionAndSideSwapTransformation = new G4RotationMatrix; reflectionAndSideSwapTransformation->rotateZ(M_PI*rad); reflectionAndSideSwapTransformation->rotateY(M_PI*rad);
  sensorhybrid_yminmax.clear();
  sensorhybrid_xpos.clear();

  // Place service hybrids first
  wedge_xpos = servicehybrid6_X/2.;
  wedge_yposmin = (wedge_Rmin + wedge_Roffset);
  wedge_yposinf = 1e-6*mm;
  {
    G4double nmodules_d = (wedge_yposmin - wedge_yposinf) / (onesensor_Y + sep_Y_module_module);
    G4int nmodules_i = std::ceil(nmodules_d);
    wedge_yposmin = ((G4double) nmodules_i) * (onesensor_Y + sep_Y_module_module) + wedge_yposinf;
  }
  while (doWedgeCloseFace){ // Loop over columns
    sensorhybrid_xpos.push_back(wedge_xpos);

    wedge_yposmax = sqrt(fabs(pow(wedge_Rmax, 2) - pow(wedge_xpos + servicehybrid6_X/2., 2)));
    G4cout << "y min/max = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;
    sensorhybrid_yminmax.push_back(std::pair<G4double, G4double>(wedge_yposmin, wedge_yposmin)); // wedge_yposmax is not exactly the location of the edge.
    std::pair<G4double, G4double>& last_sensorhybrid_yminmax = sensorhybrid_yminmax.back();

    size_t i_object=0;
    wedge_ypos = wedge_yposmin;
    while (true){
      size_t n_modules_per_side;
      if (
        (ix_service<6 && i_object==0)
        ||
        ((ix_service==3 || ix_service==4) && i_object==1)
        ||
        (ix_service==4 && i_object==2)
        ||
        ((wedge_yposmax-wedge_ypos)<servicehybrid12_Y && (wedge_yposmax-wedge_ypos)>=servicehybrid6_Y)
        ) n_modules_per_side = 3;
      else n_modules_per_side = 6;

      bool doPlaceSensorHybrid=true;

      switch (n_modules_per_side){
      case 3:
      {
        wedge_ypos += servicehybrid6_Y/2.;
        if (wedge_ypos + servicehybrid6_Y/2.>=wedge_yposmax){ doPlaceSensorHybrid = false; break; } // Breaks from the switch, not the while loop!

        // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, +(wedge_fullZ-wedge_Z)/4. - servicehybrid6_Z/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid6_X << "," << servicehybrid6_Y << "," << servicehybrid6_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 6-sensor service hybrid
        if (putServiceHybrids) BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedgeActive_Close, reflectionTransformation, relpos
        );

        wedge_ypos += servicehybrid6_Y/2.;
        last_sensorhybrid_yminmax.second += servicehybrid6_Y;
        break;
      }
      case 6:
      {
        wedge_ypos += servicehybrid12_Y/2.;
        if (wedge_ypos + servicehybrid12_Y/2.>=wedge_yposmax){ doPlaceSensorHybrid = false; break; } // Breaks from the switch, not the while loop!

        // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, +(wedge_fullZ-wedge_Z)/4. - servicehybrid12_Z/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid12_X << "," << servicehybrid12_Y << "," << servicehybrid12_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 12-sensor service hybrid
        if (putServiceHybrids) BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedgeActive_Close, reflectionTransformation, relpos
        );

        wedge_ypos += servicehybrid12_Y/2.;
        last_sensorhybrid_yminmax.second += servicehybrid12_Y;
        break;
      }
      default:
      {
        G4ExceptionDescription msg;
        msg << "n_modules_per_side = " << n_modules_per_side << " is not supported!";
        G4Exception("DetectorConstruction::DefineVolumes",
                    "MyCode0003", FatalException, msg);
        break;
      }
      }
      i_object++;

      wedge_ypos += sep_Y_servicehybrid_servicehybrid;
      if (doPlaceSensorHybrid) last_sensorhybrid_yminmax.second += sep_Y_servicehybrid_servicehybrid;
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_service++;
    last_sensorhybrid_yminmax.second -= sep_Y_servicehybrid_servicehybrid;

    wedge_xpos += (twosensor_X + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    while (
      wedge_yposmin>0.
      && (
      (wedge_xpos - servicehybrid6_X/2.)>=(wedge_Rmin + wedge_Roffset)
        ||
        (wedge_yposmin - sqrt(fabs(pow(wedge_Rmin + wedge_Roffset, 2) - pow(wedge_xpos - servicehybrid6_X/2., 2))))>=(onesensor_Y+sep_Y_module_module)
        )
      ){
      if (wedge_yposmin<(onesensor_Y + sep_Y_module_module)) break;
      G4cout << "Subtracting deltaY = " << (onesensor_Y+sep_Y_module_module) << " from y pos = " << wedge_yposmin << G4endl;
      wedge_yposmin -= (onesensor_Y + sep_Y_module_module);
      G4cout << "New y pos = " << wedge_yposmin << G4endl;
    }

    if (wedge_xpos>=wedge_Rmax) break;
  }

  // Place one or two-sensor modules next
  moduleSensorHybridConnection_yminmax_cend = sensorhybrid_yminmax.cend();
  moduleSensorHybridConnection_yminmax_left = sensorhybrid_yminmax.cbegin();
  moduleSensorHybridConnection_yminmax_right = moduleSensorHybridConnection_yminmax_left+1;
  moduleSensorHybridConnection_xpos_cend = sensorhybrid_xpos.cend();
  moduleSensorHybridConnection_xpos_left = sensorhybrid_xpos.cbegin();
  moduleSensorHybridConnection_xpos_right = moduleSensorHybridConnection_xpos_left+1;
  wedge_xpos = servicehybrid6_X + sep_X_module_servicehybrid + twosensor_X/2.;
  while (doWedgeCloseFace){ // Loop over columns
    if (moduleSensorHybridConnection_yminmax_left!=moduleSensorHybridConnection_yminmax_cend && moduleSensorHybridConnection_yminmax_right!=moduleSensorHybridConnection_yminmax_cend){
      wedge_yposmin = std::min(moduleSensorHybridConnection_yminmax_left->first, moduleSensorHybridConnection_yminmax_right->first);
      wedge_yposmax = std::max(moduleSensorHybridConnection_yminmax_left->second, moduleSensorHybridConnection_yminmax_right->second);
    }
    else if (moduleSensorHybridConnection_yminmax_left!=moduleSensorHybridConnection_yminmax_cend){
      wedge_yposmin = moduleSensorHybridConnection_yminmax_left->first;
      wedge_yposmax = moduleSensorHybridConnection_yminmax_left->second;
    }
    else if (moduleSensorHybridConnection_yminmax_right!=moduleSensorHybridConnection_yminmax_cend){
      wedge_yposmin = moduleSensorHybridConnection_yminmax_right->first;
      wedge_yposmax = moduleSensorHybridConnection_yminmax_right->second;
    }
    G4cout << "Sensor y inf/sup = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;

    bool firstColumn = (wedge_xpos<sensorhybrid_xpos.front());
    bool firstColumn_OneSensorModule = (firstColumn && sensorhybrid_xpos.front()>=onesensor_X+sep_X_module_servicehybrid);
    bool lastColumn = (wedge_xpos>sensorhybrid_xpos.back());
    bool lastColumn_OneSensorModule = (lastColumn && (wedge_Rmax-sensorhybrid_xpos.back())>=onesensor_X+sep_X_module_servicehybrid);
    bool useOneSensorModule = (firstColumn_OneSensorModule || lastColumn_OneSensorModule);
    G4double const& moduleWidthX = (useOneSensorModule ? onesensor_X : twosensor_X);
    G4double const& moduleWidthY = (useOneSensorModule ? onesensor_Y : twosensor_Y);
    G4cout << "Placing a set of " << (useOneSensorModule ? "one-" : "two-") << "sensor modules with width " << moduleWidthX << G4endl;

    size_t i_object=0;
    wedge_ypos = wedge_yposmin;
    while (true){
      G4double wedge_xpos_offset = 0; // This is an offset to calculate if placing left- or right-flanked one-sensor modules in an otherwise two-sensor module column
      bool placeLeftFlankedOneSensorModule = false;
      bool placeRightFlankedOneSensorModule = false;
      bool placeRegularModule = true;
      if (!useOneSensorModule){
        // Check if a left- or right-flanked one-sensor modules needs to be placed. If so, turn off the regular module flag.
        if (
          moduleSensorHybridConnection_yminmax_left != moduleSensorHybridConnection_yminmax_cend
          &&
          moduleSensorHybridConnection_yminmax_right != moduleSensorHybridConnection_yminmax_cend
          ){
          placeRightFlankedOneSensorModule = (
            (wedge_ypos>=moduleSensorHybridConnection_yminmax_left->second && wedge_ypos<moduleSensorHybridConnection_yminmax_right->second)
            ||
            (wedge_ypos<moduleSensorHybridConnection_yminmax_left->first && wedge_ypos>=moduleSensorHybridConnection_yminmax_right->first)
            );
          placeLeftFlankedOneSensorModule = (
            (wedge_ypos<moduleSensorHybridConnection_yminmax_left->second && wedge_ypos>=moduleSensorHybridConnection_yminmax_right->second)
            ||
            (wedge_ypos>=moduleSensorHybridConnection_yminmax_left->first && wedge_ypos<moduleSensorHybridConnection_yminmax_right->first)
            );
          placeRegularModule = !(placeRightFlankedOneSensorModule || placeLeftFlankedOneSensorModule);
        }
        // Check if the placement of a left- or right-flanekd one-sensor module makes sense.
        // Notice placeRegularModule is off if these are true.
        if (
          placeRightFlankedOneSensorModule
          &&
          std::pow(wedge_ypos, 2) + std::pow(wedge_xpos + moduleWidthX/2. - onesensor_X, 2) < pow(wedge_Rmin + wedge_Roffset, 2)
          ) placeRightFlankedOneSensorModule = false;
        if (
          placeLeftFlankedOneSensorModule
          &&
          std::pow(wedge_ypos + onesensor_Y, 2) + std::pow(wedge_xpos - moduleWidthX/2. + onesensor_X, 2) > pow(wedge_Rmax, 2)
          ) placeLeftFlankedOneSensorModule = false;
      }

      if (
        placeRegularModule
        && (
          std::pow(wedge_ypos, 2) + std::pow(wedge_xpos - moduleWidthX/2., 2) < pow(wedge_Rmin + wedge_Roffset, 2)
          ||
          std::pow(wedge_ypos + moduleWidthY, 2) + std::pow(wedge_xpos + moduleWidthX/2., 2) > pow(wedge_Rmax, 2)
          )
        ) placeRegularModule = false;

      if (placeLeftFlankedOneSensorModule) wedge_xpos_offset = (-moduleWidthX+onesensor_X)/2.;
      else if (placeRightFlankedOneSensorModule) wedge_xpos_offset = (+moduleWidthX-onesensor_X)/2.;
      if (placeLeftFlankedOneSensorModule)
        G4cout << "\t- Special left-flanked placement with an x-offset of " << wedge_xpos_offset << " relative to the column center." << G4endl;
      else if (placeRightFlankedOneSensorModule)
        G4cout << "\t- Special right-flanked placement with an x-offset of " << wedge_xpos_offset << " relative to the column center." << G4endl;
      else if (placeRegularModule)
        G4cout << "\t- No special placement; x-offset is " << wedge_xpos_offset << " relative to the column center." << G4endl;

      if ((useOneSensorModule && placeRegularModule) || placeLeftFlankedOneSensorModule || placeRightFlankedOneSensorModule){ // One-sensor modules
        wedge_ypos += onesensor_Y/2.;
        if (wedge_ypos + onesensor_Y/2.>wedge_yposmax) break; // Breaks from the switch, not the while loop!

                                                              // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos+wedge_xpos_offset, wedge_ypos, +(wedge_fullZ-wedge_Z)/4. - onesensor_Z/2.);
        G4cout << "\t- Placing one sensor modules of sizes (" << onesensor_X << "," << onesensor_Y << "," << onesensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        if (putModules) BuildOneSensorModule(
          logicWedgeActive_Close, ((lastColumn_OneSensorModule || placeLeftFlankedOneSensorModule) ? reflectionTransformation : reflectionAndSideSwapTransformation), relpos
        );

        wedge_ypos += onesensor_Y/2.;
        i_object++;
      }
      else if (!useOneSensorModule && placeRegularModule){ // Two-sensor modules
        wedge_ypos += twosensor_Y/2.;
        if (wedge_ypos + twosensor_Y/2.>wedge_yposmax) break; // Breaks from the switch, not the while loop!

                                                              // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos+wedge_xpos_offset, wedge_ypos, +(wedge_fullZ-wedge_Z)/4. - twosensor_Z/2.);
        G4cout << "\t- Placing two sensor modules of sizes (" << twosensor_X << "," << twosensor_Y << "," << twosensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        if (putModules) BuildTwoSensorModule(
          logicWedgeActive_Close, reflectionTransformation, relpos
        );

        wedge_ypos += twosensor_Y/2.;
        i_object++;
      }
      else wedge_ypos += (useOneSensorModule ? onesensor_Y : twosensor_Y);

      wedge_ypos += sep_Y_module_module;
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_module++;

    wedge_xpos += (moduleWidthX/2. + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    G4double const& next_moduleWidthX = (wedge_xpos<sensorhybrid_xpos.back() ? twosensor_X : onesensor_X);
    wedge_xpos += next_moduleWidthX/2.;

    if (wedge_xpos+next_moduleWidthX/2.>=wedge_Rmax) break;
    if (moduleSensorHybridConnection_yminmax_right==moduleSensorHybridConnection_yminmax_cend) break;

    moduleSensorHybridConnection_yminmax_left = moduleSensorHybridConnection_yminmax_right;
    moduleSensorHybridConnection_yminmax_right++;
    moduleSensorHybridConnection_xpos_left = moduleSensorHybridConnection_xpos_right;
    moduleSensorHybridConnection_xpos_right++;
  }


  /**********************************/
  /**********************************/
  /* Construct the wedge components */
  /**********************************/
  /**********************************/
  if (putWedgeComponents) BuildWedgeComponents(logicWedgePassive, logicWedgeActive_Far, coolingpipes_xpos_ymin);
  if (doEndCapSupport) BuildEndCapSupportComponents(logicEndcap);

  // Place the passive material
  detname = det_wedge+"_PassiveMaterial";
  new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logicWedgePassive,
    detname.c_str(),
    logicWedge,
    false,
    0,
    fCheckOverlaps
  );
  // Place the active material
  detname = det_wedge+"_ActiveMaterial";
  new G4PVPlacement(
    nullptr,
    G4ThreeVector(0, 0, -(wedge_fullZ+wedge_Z)/4.),
    logicWedgeActive_Close,
    (detname+"_Close").c_str(),
    logicWedge,
    false,
    0,
    fCheckOverlaps
  );
  new G4PVPlacement(
    nullptr,
    G4ThreeVector(0, 0, +(wedge_fullZ+wedge_Z)/4.),
    logicWedgeActive_Far,
    (detname+"_Far").c_str(),
    logicWedge,
    false,
    0,
    fCheckOverlaps
  );
  // Place the disks
  detname = det_disk;
  if (doCloseDisk){
    G4RotationMatrix* rotateDisk = new G4RotationMatrix; rotateDisk->rotateY(M_PI*rad);
    new G4PVPlacement(
      rotateDisk,
      G4ThreeVector(0, 0, -(endcap_Z-diskBox_Z)/2.),
      logicDiskBox,
      (detname+"_Close").c_str(),
      logicEndcap,
      false,
      0,
      fCheckOverlaps
    );
  }
  if (doFarDisk){
    new G4PVPlacement(
      nullptr,
      G4ThreeVector(0, 0, (endcap_Z-diskBox_Z)/2.),
      logicDiskBox,
      (detname+"_Far").c_str(),
      logicEndcap,
      false,
      0,
      fCheckOverlaps
    );
  }

  // Place the end caps
  // IMPORTANT NOTE: reflections have to come at the very last step; otherwise nothing gets reflected at all!
  detname = det_endcap;
  if (!endcapsAreReflected){
    if (doBackEndcap){
      G4RotationMatrix* rotateETL = new G4RotationMatrix; rotateETL->rotateY(M_PI*rad);
      new G4PVPlacement(
        rotateETL,
        G4ThreeVector(0, 0, -(diskPair_Z + IP_dZ)/2.),
        logicEndcap,
        (detname+"_Back").c_str(),
        logicMTD,
        false,
        0,
        fCheckOverlaps
      );
    }
  }
  else{
    if (doBackEndcap){
      G4Translate3D translation(0, 0, -(diskPair_Z + IP_dZ)/2.);
      G4RotationMatrix* rot3D = new G4RotationMatrix();
      G4Transform3D rotation = G4Rotate3D(*rot3D);
      G4ReflectZ3D reflection;
      G4Transform3D tr3D = translation*rotation*reflection;
      G4ReflectionFactory::Instance()->Place(
        tr3D,
        (detname+"_Back").c_str(),
        logicEndcap,
        logicMTD,
        false,
        0,
        fCheckOverlaps
      );
    }
  }
  if (doFrontEndcap) new G4PVPlacement(
    nullptr,
    G4ThreeVector(0, 0, (diskPair_Z + IP_dZ)/2.),
    logicEndcap,
    (detname+"_Front").c_str(),
    logicMTD,
    false,
    0,
    fCheckOverlaps
  );

  fScoringVolume = logicWedge;

  //always return the physical World
  return physWorld;
}
