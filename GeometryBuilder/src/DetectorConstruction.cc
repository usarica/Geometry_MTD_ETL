#include <cassert>
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
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"


using namespace CLHEP;
using namespace std;


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
  G4double a;
  G4double z;
  G4double density;

  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* CO2 = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  //nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

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

  // When calling FindOrBuildElement, do not use 'G4_' at the beginning of the element name
  G4Element* Aluminum = nist->FindOrBuildElement("Al");
  G4Element* Nitrogen = nist->FindOrBuildElement("N");
  G4Element* Carbon = nist->FindOrBuildElement("C");
  G4Element* Hydrogen = nist->FindOrBuildElement("H");
  G4Element* Silicon = nist->FindOrBuildElement("Si");
  G4Element* Tin = nist->FindOrBuildElement("Sn");
  G4Element* Silver = nist->FindOrBuildElement("Ag");
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

  G4Material* FR4 = new G4Material("FR4", 1.85*g/cm3, 1, kStateSolid, NTP_Temperature, STP_Pressure); // For circuit boards, density approximate?
  FR4->AddElement(Carbon, (G4int) 1);

  G4Material* PCB = new G4Material("PCB", 1.85*g/cm3, 1, kStateSolid, NTP_Temperature, STP_Pressure); // For circuit boards, density approximate?
  PCB->AddElement(Carbon, (G4int) 1);

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

  G4Material* SnAg = new G4Material("SnAg", (7.31*0.5 + 10.49*0.5)*g/cm3, 2, kStateSolid, NTP_Temperature, STP_Pressure); // FIXME: Composition should be per mol, and density needs recalc.
  SnAg->AddElement(Tin, (G4double) 0.5);
  SnAg->AddElement(Silver, (G4double) 0.5);

  G4Material* SensorBump = new G4Material("SensorBump", (SnAg->GetDensity()*0.1 + Air->GetDensity()*0.9), 2, kStateSolid, NTP_Temperature, STP_Pressure); // FIXME: Composition should be per mol, and density needs recalc
  SensorBump->AddMaterial(SnAg, (G4double) 0.1);
  SensorBump->AddMaterial(Air, (G4double) 0.9);

  G4Material* ServiceConnector = new G4Material("ServiceConnector", (8.96*0.1 + Air->GetDensity()*0.9), 2, kStateSolid, NTP_Temperature, STP_Pressure); // FIXME: Composition should be per mol, and density needs recalc
  ServiceConnector->AddElement(Copper, (G4double) 0.1);
  ServiceConnector->AddMaterial(Air, (G4double) 0.9);

  G4cout << "DetectorConstruction::DefineMaterials: AlN built!" << G4endl;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


void DetectorConstruction::BuildOneSensorModule(
  bool rightFlank,
  G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos,
  G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV
){
  using namespace ETLDetectorDimensions;

  // See Fig. 3.56 in the MTD TDR for the module diagrams
  // Also Fig. 3.73 for the z positions.
  // Note that the base aluminum plate on which the service is loaded, or the epoxy underneath are not included since they are only on one side.
  string const detbase = "ETLOneSensorModule";
  string detname;

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
  G4Material* lgadMat = G4Material::GetMaterial("Mat_LGAD");
  G4VisAttributes* lgadVisAttr = new G4VisAttributes(G4Colour::Yellow()); lgadVisAttr->SetVisibility(true);

  // Solder bumps
  detname = detbase + "_SolderBumps";
  G4double bumpsSize_X = getDimension(detname+"_X");
  G4double bumpsSize_Y = getDimension(detname+"_Y");
  G4double bumpsSize_Z = getDimension(detname+"_Z");
  G4Material* bumpsMat = G4Material::GetMaterial("SensorBump");
  G4VisAttributes* bumpsVisAttr = new G4VisAttributes(G4Colour::Magenta()); bumpsVisAttr->SetVisibility(true);

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  G4double epoxySize_X = getDimension(detname+"_X");
  G4double epoxySize_Y = getDimension(detname+"_Y");
  G4double epoxySize_Z = getDimension(detname+"_Z");
  G4Material* epoxyMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* epoxyVisAttr = new G4VisAttributes(G4Colour::Brown()); epoxyVisAttr->SetVisibility(true);

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  G4double coverplateSize_X = getDimension(detname+"_X");
  G4double coverplateSize_Y = getDimension(detname+"_Y");
  G4double coverplateSize_Z = getDimension(detname+"_Z");
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
  moduleBox = new G4Box(
    detname.c_str(),
    moduleSize_X/2., moduleSize_Y/2., moduleSize_Z/2.
  );
  moduleLogical = new G4LogicalVolume(
    moduleBox, 
    moduleMat,
    detname.c_str()
  );
  moduleLogical->SetVisAttributes(moduleVisAttr);
  modulePV = new G4PVPlacement(
    0, relativePos,
    moduleLogical,
    detname.c_str(),
    motherLogical,
    false, 0, fCheckOverlaps
  );

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
  G4PVPlacement* basefilmPV = new G4PVPlacement(
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
  G4PVPlacement* baseplatePV = new G4PVPlacement(
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
  G4PVPlacement* lairdfilmPV = new G4PVPlacement(
    0, G4ThreeVector((lairdfilmSize_X-baseplateSize_X)/2.*(rightFlank ? 1. : -1.), 0, basefilmSize_Z+baseplateSize_Z+(lairdfilmSize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* etroc1PV = new G4PVPlacement(
    0, G4ThreeVector((etrocSize_X-baseplateSize_X)/2.*(rightFlank ? 1. : -1.), -etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc2PV = new G4PVPlacement(
    0, G4ThreeVector((etrocSize_X-baseplateSize_X)/2.*(rightFlank ? 1. : -1.), +etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* bumpsPV = new G4PVPlacement(
    0, G4ThreeVector((bumpsSize_X-baseplateSize_X)/2.*(rightFlank ? 1. : -1.), 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    detname.c_str(),
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
  G4PVPlacement* lgadPV = new G4PVPlacement(
    0, G4ThreeVector((lgadSize_X-baseplateSize_X)/2.*(rightFlank ? 1. : -1.), 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* epoxyPV = new G4PVPlacement(
    0, G4ThreeVector((epoxySize_X-baseplateSize_X)/2.*(rightFlank ? 1. : -1.), 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* coverplatePV = new G4PVPlacement(
    0, G4ThreeVector((coverplateSize_X-baseplateSize_X)/2.*(rightFlank ? 1. : -1.), 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical,
    detname.c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
}

void DetectorConstruction::BuildTwoSensorModule(
  G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos,
  G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV
){
  using namespace ETLDetectorDimensions;

  // See Fig. 3.56 in the MTD TDR for the module diagrams
  // Also Fig. 3.73 for the z positions.
  // Note that the base aluminum plate on which the service is loaded, or the epoxy underneath are not included since they are only on one side.
  string const detbase = "ETLTwoSensorModule";
  string detname;

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
  G4Material* lgadMat = G4Material::GetMaterial("Mat_LGAD");
  G4VisAttributes* lgadVisAttr = new G4VisAttributes(G4Colour::Yellow()); lgadVisAttr->SetVisibility(true);

  // Solder bumps
  detname = detbase + "_SolderBumps";
  G4double bumpsSize_X = getDimension(detname+"_X");
  G4double bumpsSize_Y = getDimension(detname+"_Y");
  G4double bumpsSize_Z = getDimension(detname+"_Z");
  G4Material* bumpsMat = G4Material::GetMaterial("SensorBump");
  G4VisAttributes* bumpsVisAttr = new G4VisAttributes(G4Colour::Magenta()); bumpsVisAttr->SetVisibility(true);

  // Epoxy betwen sensor and AlN cover
  detname = detbase + "_EpoxyCover";
  G4double epoxySize_X = getDimension(detname+"_X");
  G4double epoxySize_Y = getDimension(detname+"_Y");
  G4double epoxySize_Z = getDimension(detname+"_Z");
  G4Material* epoxyMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* epoxyVisAttr = new G4VisAttributes(G4Colour::Brown()); epoxyVisAttr->SetVisibility(true);

  // AlN sensor cover
  detname = detbase + "_CoverPlate";
  G4double coverplateSize_X = getDimension(detname+"_X");
  G4double coverplateSize_Y = getDimension(detname+"_Y");
  G4double coverplateSize_Z = getDimension(detname+"_Z");
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
  moduleBox = new G4Box(
    detname.c_str(),
    moduleSize_X/2., moduleSize_Y/2., moduleSize_Z/2.
  );
  moduleLogical = new G4LogicalVolume(
    moduleBox,
    moduleMat,
    detname.c_str()
  );
  moduleLogical->SetVisAttributes(moduleVisAttr);
  modulePV = new G4PVPlacement(
    0, relativePos,
    moduleLogical,
    detname.c_str(),
    motherLogical,
    false, 0, fCheckOverlaps
  );

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
  G4PVPlacement* basefilmPV = new G4PVPlacement(
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
  G4PVPlacement* baseplatePV = new G4PVPlacement(
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
  G4PVPlacement* lairdfilmPV = new G4PVPlacement(
    0, G4ThreeVector((lairdfilmSize_X-baseplateSize_X)/2., 0, basefilmSize_Z+baseplateSize_Z+(lairdfilmSize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* etroc1PV = new G4PVPlacement(
    0, G4ThreeVector(-etrocSize_X/2., -etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc2PV = new G4PVPlacement(
    0, G4ThreeVector(-etrocSize_X/2., +etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc3PV = new G4PVPlacement(
    0, G4ThreeVector(+etrocSize_X/2., -etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical,
    (detname+"3").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc4PV = new G4PVPlacement(
    0, G4ThreeVector(+etrocSize_X/2., +etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* bumps1PV = new G4PVPlacement(
    0, G4ThreeVector(-bumpsSize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* bumps2PV = new G4PVPlacement(
    0, G4ThreeVector(+bumpsSize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* lgad1PV = new G4PVPlacement(
    0, G4ThreeVector(-lgadSize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
    lgadLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* lgad2PV = new G4PVPlacement(
    0, G4ThreeVector(+lgadSize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* epoxy1PV = new G4PVPlacement(
    0, G4ThreeVector(-epoxySize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
    epoxyLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* epoxy2PV = new G4PVPlacement(
    0, G4ThreeVector(+epoxySize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
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
  G4PVPlacement* coverplate1PV = new G4PVPlacement(
    0, G4ThreeVector(-coverplateSize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical,
    (detname+"1").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* coverplate2PV = new G4PVPlacement(
    0, G4ThreeVector(+coverplateSize_X/2., 0, basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical,
    (detname+"2").c_str(),
    moduleLogical,
    false, 0, fCheckOverlaps
  );
}

void DetectorConstruction::BuildSensorServiceHybrid(
  int const& nSensorsPerSide, // 6 or 3
  G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos,
  G4Box*& serviceBox, G4LogicalVolume*& serviceLogical, G4PVPlacement*& servicePV
){
  using namespace ETLDetectorDimensions;

  string const detbase = "ETL" + std::to_string(2*nSensorsPerSide) + "SensorServiceHybrid";
  string detname;

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


  // Build the module
  detname = detbase;
  serviceBox = new G4Box(
    detname.c_str(),
    servicehybridSize_X/2., servicehybridSize_Y/2., servicehybridSize_Z/2.
  );
  serviceLogical = new G4LogicalVolume(
    serviceBox,
    servicehybridMat,
    detname.c_str()
  );
  serviceLogical->SetVisAttributes(servicehybridVisAttr);
  servicePV = new G4PVPlacement(
    0, relativePos,
    serviceLogical,
    detname.c_str(),
    motherLogical,
    false, 0, fCheckOverlaps
  );


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
  G4PVPlacement* thermalpadPV = new G4PVPlacement(
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
  G4PVPlacement* readoutboardPV = new G4PVPlacement(
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
  G4PVPlacement* connectorPV = new G4PVPlacement(
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
  G4PVPlacement* powerboardPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, thermalpadSize_Z+readoutboardSize_Z+connectorSize_Z+(powerboardSize_Z-servicehybridSize_Z)/2.),
    powerboardLogical,
    detname.c_str(),
    serviceLogical,
    false, 0, fCheckOverlaps
  );
}



G4VPhysicalVolume* DetectorConstruction::DefineVolumes(){
  using namespace ETLDetectorDimensions;

  // Setup materials
  G4Material* world_mat = G4Material::GetMaterial("Vacuum");
  G4Material* envelope_mat = G4Material::GetMaterial("Vacuum");
  G4Material* wedge_mat = G4Material::GetMaterial("AlN");
  if (!world_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve world_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }
  if (!envelope_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve envelope_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }
  if (!wedge_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve wedge_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }

  // Size parameters
  string const det_offset = "ETLOffset";
  string const det_wedge = "ETLWedge";
  string const det_onesensor = "ETLOneSensorModule";
  string const det_twosensor = "ETLTwoSensorModule";
  string const det_servicehybrid6 = "ETL6SensorServiceHybrid";
  string const det_servicehybrid12 = "ETL12SensorServiceHybrid";

  string detname;

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
  G4double wedge_fullZ = wedge_Z + 2.*std::max(std::max(servicehybrid12_Z, servicehybrid6_Z), std::max(onesensor_Z, twosensor_Z));

  // External offsets
  detname = det_offset + "_Module_SensorServiceHybrid_dX";
  G4double sep_X_module_servicehybrid = getDimension(detname);
  detname = det_offset + "_Module_Module_dY";
  G4double sep_Y_module_module = getDimension(detname);
  detname = det_offset + "_SensorServiceHybrid_SensorServiceHybrid_dY";
  G4double sep_Y_servicehybrid_servicehybrid = getDimension(detname);

  // World parameters
  G4double world_sizeX = wedge_Rmax*2.;
  G4double world_sizeY = wedge_Rmax*2.;
  G4double world_sizeZ = wedge_fullZ;

  // Generic envelope parameters
  G4double env_sizeX = world_sizeX;
  G4double env_sizeY = world_sizeY;
  G4double env_sizeZ = world_sizeZ;


  // World
  G4Box* solidWorld = new G4Box(
    "World",
    0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ
  );
  G4LogicalVolume* logicWorld = new G4LogicalVolume(
    solidWorld,          //its solid
    world_mat,           //its material
    "World"              //its name
  );
  G4VisAttributes* worldVisAttr = new G4VisAttributes(G4Colour::Black()); worldVisAttr->SetVisibility(false);
  logicWorld->SetVisAttributes(worldVisAttr);
  G4VPhysicalVolume* physWorld = new G4PVPlacement(
    0,                     //no rotation
    G4ThreeVector(),       //at (0,0,0)
    logicWorld,            //its logical volume
    "World",               //its name
    0,                     //its mother  volume
    false,                 //no boolean operation
    0,                     //copy number
    fCheckOverlaps         //overlaps checking
  );

  // Envelope
  G4Box* solidEnvelope = new G4Box(
    "Envelope",
    0.5*env_sizeX, 0.5*env_sizeY, 0.5*env_sizeZ
  );
  G4LogicalVolume* logicEnvelope = new G4LogicalVolume(
    solidEnvelope,          //its solid
    envelope_mat,           //its material
    "Envelope"              //its name
  );
  G4VisAttributes* envVisAttr = new G4VisAttributes(G4Colour::Black()); envVisAttr->SetVisibility(false);
  logicEnvelope->SetVisAttributes(envVisAttr);
  new G4PVPlacement(
    0,                       //no rotation
    G4ThreeVector(),         //at (0,0,0)
    logicEnvelope,           //its logical volume
    "Envelope",              //its name
    logicWorld,              //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    fCheckOverlaps           //overlaps checking
  );

  // Wedge
  detname = det_wedge;
  G4Tubs* solidWedge = new G4Tubs(
    detname.c_str(),
    wedge_Rmin, wedge_Rmax, wedge_fullZ/2., 0, 90.*degree
  );
  G4LogicalVolume* logicWedge = new G4LogicalVolume(
    solidWedge,
    wedge_mat,
    detname.c_str()
  );
  G4VisAttributes* wedgeVisAttr = new G4VisAttributes(G4Colour::Gray()); wedgeVisAttr->SetVisibility(true);
  logicWedge->SetVisAttributes(wedgeVisAttr);
  new G4PVPlacement(
    0,                       //no rotation
    G4ThreeVector(),         //at (0,0,0)
    logicWedge,                //its logical volume
    detname.c_str(),              //its name
    logicEnvelope,              //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    fCheckOverlaps           //overlaps checking
  );


  /*****************************************/
  /*****************************************/
  /* Construct the front face of the wedge */
  /*****************************************/
  /*****************************************/
  size_t ix_module, ix_service;
  G4double wedge_xpos;
  G4double wedge_ypos;
  G4double wedge_yposmin;
  G4double wedge_yposmax;
  G4double wedge_Roffset;
  std::vector<std::pair<G4double, G4double>> sensorhybrid_yminmax;
  std::pair<G4double, G4double>::const_iterator moduleSensorHybridConnection_left;
  std::pair<G4double, G4double>::const_iterator moduleSensorHybridConnection_right;

  wedge_Roffset=0;
  ix_module=ix_service=0;

  // Place service hybrids first
  wedge_xpos = onesensor_X + sep_X_module_servicehybrid + servicehybrid6_X/2.;
  G4double first_sensorhybrid_x = wedge_xpos;
  wedge_yposmin = (wedge_Rmin + wedge_Roffset);
  while (true){ // Loop over columns
    wedge_yposmax = sqrt(fabs(pow(wedge_Rmax, 2) - pow(wedge_xpos + servicehybrid6_X/2., 2)));
    G4cout << "y min/max = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;
    sensorhybrid_yminmax.push_back(std::pair<G4double, G4double>(wedge_yposmin, wedge_yposmax));

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

      G4Box* servicehybridBox = nullptr;
      G4LogicalVolume* servicehybridLogical = nullptr;
      G4PVPlacement* servicehybridPV = nullptr;

      switch (n_modules_per_side){
      case 3:
      {
        wedge_ypos += servicehybrid6_Y/2.;
        if (wedge_ypos + servicehybrid6_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

        // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+servicehybrid6_Z)/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid6_X << "," << servicehybrid6_Y << "," << servicehybrid6_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 6-sensor service hybrid
        BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedge, relpos,
          servicehybridBox, servicehybridLogical, servicehybridPV
        );

        wedge_ypos += servicehybrid6_Y/2.;
        break;
      }
      case 6:
      {
        wedge_ypos += servicehybrid12_Y/2.;
        if (wedge_ypos + servicehybrid12_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

        // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+servicehybrid12_Z)/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid12_X << "," << servicehybrid12_Y << "," << servicehybrid12_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 12-sensor service hybrid
        BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedge, relpos,
          servicehybridBox, servicehybridLogical, servicehybridPV
        );

        wedge_ypos += servicehybrid12_Y/2.;
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
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_service++;

    wedge_xpos += (twosensor_X + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    while (
      wedge_yposmin>0.
      && (
        (wedge_xpos - servicehybrid6_X/2.)>=(wedge_Rmin + wedge_Roffset)
        ||
        (wedge_yposmin - sqrt(fabs(pow(wedge_Rmin + wedge_Roffset, 2) - pow(wedge_xpos - servicehybrid6_X/2., 2))))>=(onesensor_Y+sep_Y_module_module)
        )
      ){
      if (wedge_yposmin<(onesensor_Y+sep_Y_module_module)) break;
      G4cout << "Subtracting deltaY = " << (onesensor_Y+sep_Y_module_module) << " from y pos = " << wedge_yposmin << G4endl;
      wedge_yposmin -= (onesensor_Y+sep_Y_module_module);
    }

    if (wedge_xpos>=wedge_Rmax) break;
  }
  G4double last_sensorhybrid_x = wedge_xpos;

  // Place one or two-sensor modules next
  wedge_xpos = onesensor_X/2.;
  moduleSensorHybridConnection_left = sensorhybrid_yminmax.cend();
  moduleSensorHybridConnection_right = sensorhybrid_yminmax.cbegin();
  wedge_yposmin = *moduleSensorHybridConnection_right;
  while (true){ // Loop over columns
    bool firstColumn = (wedge_xpos<first_sensorhybrid_x);
    bool lastColumn = (wedge_xpos>last_sensorhybrid_x);

    bool useOneSensorModule = (firstColumn || lastColumn);
    G4double const& moduleWidthX = (useOneSensorModule ? onesensor_X : twosensor_X);
    G4cout << "Placing a set of " << (useOneSensorModule ? "one-" : "two-") << "sensor modules with width " << moduleWidthX << G4endl;

    wedge_yposmax = sqrt(fabs(pow(wedge_Rmax, 2) - pow(wedge_xpos + moduleWidthX/2., 2)));
    G4cout << "Sensor y min/max = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;

    size_t i_object=0;
    wedge_ypos = wedge_yposmin;
    while (true){
      G4Box* moduleBox = nullptr;
      G4LogicalVolume* moduleLogical = nullptr;
      G4PVPlacement* modulePV = nullptr;

      if (!useOneSensorModule){ // Two-sensor modules
        wedge_ypos += twosensor_Y/2.;
        if (wedge_ypos + twosensor_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

        // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+twosensor_Z)/2.);
        G4cout << "\t- Placing two sensor modules of sizes (" << twosensor_X << "," << twosensor_Y << "," << twosensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        BuildTwoSensorModule(
          logicWedge, relpos,
          moduleBox, moduleLogical, modulePV
        );

        wedge_ypos += twosensor_Y/2.;
      }
      else{ // One-sensor modules
        wedge_ypos += onesensor_Y/2.;
        if (wedge_ypos + onesensor_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

        // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+onesensor_Z)/2.);
        G4cout << "\t- Placing one sensor modules of sizes (" << onesensor_X << "," << onesensor_Y << "," << onesensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        BuildOneSensorModule(
          firstColumn,
          logicWedge, relpos,
          moduleBox, moduleLogical, modulePV
        );

        wedge_ypos += onesensor_Y/2.;
      }
      i_object++;

      wedge_ypos += sep_Y_module_module;
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_module++;

    wedge_xpos += (moduleWidthX/2. + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    while (
      wedge_yposmin>0.
      && (
      wedge_xpos>=(wedge_Rmin + wedge_Roffset)
        ||
        (wedge_yposmin - sqrt(fabs(pow(wedge_Rmin + wedge_Roffset, 2) - pow(wedge_xpos, 2))))>=(onesensor_Y+sep_Y_module_module)
        )
      ){
      if (wedge_yposmin<(onesensor_Y+sep_Y_module_module)) break;
      //G4cout << "Subtracting deltaY = " << (onesensor_Y+sep_Y_module_module) << " from y pos = " << wedge_yposmin << G4endl;
      wedge_yposmin -= (onesensor_Y+sep_Y_module_module);
    }
    G4double const& next_moduleWidthX = (wedge_xpos<last_sensorhybrid_x ? twosensor_X : onesensor_X);
    wedge_xpos += next_moduleWidthX/2.;

    if (wedge_xpos+next_moduleWidthX/2.>=wedge_Rmax) break;
  }


  /****************************************/
  /****************************************/
  /* Construct the back face of the wedge */
  /****************************************/
  /****************************************/
  wedge_Roffset=0;
  ix_module=ix_service=0;
  wedge_xpos=0;
  sensorhybrid_yminmax.clear();

  // Place service hybrids first
  wedge_xpos = servicehybrid6_X/2.;
  G4double first_sensorhybrid_x = wedge_xpos;
  wedge_yposmin = (wedge_Rmin + wedge_Roffset);
  while (true){ // Loop over columns
    wedge_yposmax = sqrt(fabs(pow(wedge_Rmax, 2) - pow(wedge_xpos + servicehybrid6_X/2., 2)));
    G4cout << "y min/max = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;
    sensorhybrid_yminmax.push_back(std::pair<G4double, G4double>(wedge_yposmin, wedge_yposmax));

    size_t i_object=0;
    wedge_ypos = wedge_yposmin;
    while (true){
      size_t n_modules_per_side;
      if (
        (ix_service<6 && i_object==0)
        ||
        (ix_service==3 && i_object<2)
        ||
        (ix_service==4 && i_object<3)
        ||
        ((wedge_yposmax-wedge_ypos)<servicehybrid12_Y && (wedge_yposmax-wedge_ypos)>=servicehybrid6_Y)
        ) n_modules_per_side = 3;
      else n_modules_per_side = 6;

      G4Box* servicehybridBox = nullptr;
      G4LogicalVolume* servicehybridLogical = nullptr;
      G4PVPlacement* servicehybridPV = nullptr;

      switch (n_modules_per_side){
      case 3:
      {
        wedge_ypos += servicehybrid6_Y/2.;
        if (wedge_ypos + servicehybrid6_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

                                                                    // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+servicehybrid6_Z)/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid6_X << "," << servicehybrid6_Y << "," << servicehybrid6_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 6-sensor service hybrid
        BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedge, relpos,
          servicehybridBox, servicehybridLogical, servicehybridPV
        );

        wedge_ypos += servicehybrid6_Y/2.;
        break;
      }
      case 6:
      {
        wedge_ypos += servicehybrid12_Y/2.;
        if (wedge_ypos + servicehybrid12_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

                                                                     // Position of service hybrid center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+servicehybrid12_Z)/2.);
        G4cout << "\t- Placing service hybrid of sizes (" << servicehybrid12_X << "," << servicehybrid12_Y << "," << servicehybrid12_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the 12-sensor service hybrid
        BuildSensorServiceHybrid(
          n_modules_per_side,
          logicWedge, relpos,
          servicehybridBox, servicehybridLogical, servicehybridPV
        );

        wedge_ypos += servicehybrid12_Y/2.;
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
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_service++;

    wedge_xpos += (twosensor_X + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    while (
      wedge_yposmin>0.
      && (
      (wedge_xpos - servicehybrid6_X/2.)>=(wedge_Rmin + wedge_Roffset)
        ||
        (wedge_yposmin - sqrt(fabs(pow(wedge_Rmin + wedge_Roffset, 2) - pow(wedge_xpos - servicehybrid6_X/2., 2))))>=(onesensor_Y+sep_Y_module_module)
        )
      ){
      if (wedge_yposmin<(onesensor_Y+sep_Y_module_module)) break;
      G4cout << "Subtracting deltaY = " << (onesensor_Y+sep_Y_module_module) << " from y pos = " << wedge_yposmin << G4endl;
      wedge_yposmin -= (onesensor_Y+sep_Y_module_module);
    }

    if (wedge_xpos>=wedge_Rmax) break;
  }
  G4double last_sensorhybrid_x = wedge_xpos;

  // Place one or two-sensor modules next
  wedge_xpos = onesensor_X/2.;
  wedge_yposmin = (wedge_Rmin + wedge_Roffset);
  moduleSensorHybridConnection_left = sensorhybrid_yminmax.cend();
  moduleSensorHybridConnection_right = sensorhybrid_yminmax.cbegin();
  while (true){ // Loop over columns
    bool firstColumn = (wedge_xpos<first_sensorhybrid_x);
    bool lastColumn = (wedge_xpos>last_sensorhybrid_x);
    bool useOneSensorModule = (firstColumn || lastColumn);
    G4double const& moduleWidthX = (useOneSensorModule ? onesensor_X : twosensor_X);
    G4cout << "Placing a set of " << (useOneSensorModule ? "one-" : "two-") << "sensor modules with width " << moduleWidthX << G4endl;

    wedge_yposmax = sqrt(fabs(pow(wedge_Rmax, 2) - pow(wedge_xpos + moduleWidthX/2., 2)));
    G4cout << "Sensor y min/max = " << wedge_yposmin << " / " << wedge_yposmax << G4endl;

    size_t i_object=0;
    wedge_ypos = wedge_yposmin;
    while (true){
      G4Box* moduleBox = nullptr;
      G4LogicalVolume* moduleLogical = nullptr;
      G4PVPlacement* modulePV = nullptr;

      if (!useOneSensorModule){ // Two-sensor modules
        wedge_ypos += twosensor_Y/2.;
        if (wedge_ypos + twosensor_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

                                                               // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+twosensor_Z)/2.);
        G4cout << "\t- Placing two sensor modules of sizes (" << twosensor_X << "," << twosensor_Y << "," << twosensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        BuildTwoSensorModule(
          logicWedge, relpos,
          moduleBox, moduleLogical, modulePV
        );

        wedge_ypos += twosensor_Y/2.;
      }
      else{ // One-sensor modules
        wedge_ypos += onesensor_Y/2.;
        if (wedge_ypos + onesensor_Y/2.>=wedge_yposmax) break; // Breaks from the switch, not the while loop!

                                                               // Position of module center relative to the wedge center
        G4ThreeVector relpos(wedge_xpos, wedge_ypos, (wedge_Z+onesensor_Z)/2.);
        G4cout << "\t- Placing one sensor modules of sizes (" << onesensor_X << "," << onesensor_Y << "," << onesensor_Z << ") at position (" << relpos.x() << "," << relpos.y() << "," << relpos.z() << ")" << G4endl;
        // Construct the module
        BuildOneSensorModule(
          firstColumn,
          logicWedge, relpos,
          moduleBox, moduleLogical, modulePV
        );

        wedge_ypos += onesensor_Y/2.;
      }
      i_object++;

      wedge_ypos += sep_Y_module_module;
      if (wedge_ypos>=wedge_yposmax) break;
    }
    ix_module++;

    wedge_xpos += (moduleWidthX/2. + sep_X_module_servicehybrid*2. + servicehybrid6_X);
    while (
      wedge_yposmin>0.
      && (
        wedge_xpos>=(wedge_Rmin + wedge_Roffset)
        ||
        (wedge_yposmin - sqrt(fabs(pow(wedge_Rmin + wedge_Roffset, 2) - pow(wedge_xpos, 2))))>=(onesensor_Y+sep_Y_module_module)
        )
      ){
      if (wedge_yposmin<(onesensor_Y+sep_Y_module_module)) break;
      //G4cout << "Subtracting deltaY = " << (onesensor_Y+sep_Y_module_module) << " from y pos = " << wedge_yposmin << G4endl;
      wedge_yposmin -= (onesensor_Y+sep_Y_module_module);
    }
    G4double const& next_moduleWidthX = (wedge_xpos<last_sensorhybrid_x ? twosensor_X : onesensor_X);
    wedge_xpos += next_moduleWidthX/2.;

    if (wedge_xpos+next_moduleWidthX/2.>=wedge_Rmax) break;
  }





  fScoringVolume = logicWedge;

  //always return the physical World
  return physWorld;
}
