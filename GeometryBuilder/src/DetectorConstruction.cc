#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"


using namespace CLHEP;


DetectorConstruction::DetectorConstruction() :
  G4VUserDetectorConstruction(),
  fCheckOverlaps(true),
  fScoringVolume(nullptr)
{}
DetectorConstruction::~DetectorConstruction(){}


G4VPhysicalVolume* DetectorConstruction::Construct(){
  G4cout << "Begin DetectorConstruction::Construct" << G4endl;
  DefineMaterials();
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
  ServiceConnector->AddMaterial(Copper, (G4double) 0.1);
  ServiceConnector->AddMaterial(Air, (G4double) 0.9);

  G4cout << "DetectorConstruction::DefineMaterials: AlN built!" << G4endl;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


void DetectorConstruction::BuildOneSensorModule(
  G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos,
  G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV
){
  // See Fig. 3.56 in the MTD TDR for the module diagrams
  // Also Fig. 3.73 for the z positions.
  // Note that the base aluminum plate on which the service is loaded, or the epoxy underneath are not included since they are only on one side.
  // AlN
  G4double baseplateSize_X = 43.1*mm;
  G4double baseplateSize_Y = 28.25*mm;
  G4double baseplateSize_Z = 0.79*mm;
  G4Material* baseplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* baseplateVisAttr = new G4VisAttributes(G4Colour::Gray()); baseplateVisAttr->SetVisibility(true);

  // Thermal pad
  G4double basefilmSize_X = baseplateSize_X;
  G4double basefilmSize_Y = baseplateSize_Y;
  G4double basefilmSize_Z = 0.25*mm;
  G4Material* basefilmMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* basefilmVisAttr = new G4VisAttributes(G4Colour::Brown()); basefilmVisAttr->SetVisibility(true);

  // ETROC
  G4double etrocSize_X = 20.8*mm;
  G4double etrocSize_Y = 22.3*mm;
  G4double etrocSize_Z = 0.25*mm;
  G4Material* etrocMat = G4Material::GetMaterial("Si");
  G4VisAttributes* etrocVisAttr = new G4VisAttributes(G4Colour::Green()); etrocVisAttr->SetVisibility(true);

  // Laird film
  G4double lairdfilmSize_X = etrocSize_X*2.;
  G4double lairdfilmSize_Y = etrocSize_Y;
  G4double lairdfilmSize_Z = 0.08*mm;
  G4Material* lairdfilmMat = G4Material::GetMaterial("Laird");
  G4VisAttributes* lairdfilmVisAttr = new G4VisAttributes(G4Colour::Brown()); lairdfilmVisAttr->SetVisibility(true);

  // LGAD sensor
  G4double lgadSize_X = 42.0*mm;
  G4double lgadSize_Y = 21.2*mm;
  G4double lgadSize_Z = 0.3*mm;
  G4Material* lgadMat = G4Material::GetMaterial("Si");
  G4VisAttributes* lgadVisAttr = new G4VisAttributes(G4Colour::Yellow()); lgadVisAttr->SetVisibility(true);

  // Solder bumps
  G4double bumpsSize_X = etrocSize_X*2.;
  G4double bumpsSize_Y = lgadSize_Y;
  G4double bumpsSize_Z = 0.03*mm;
  G4Material* bumpsMat = G4Material::GetMaterial("SensorBump");
  G4VisAttributes* bumpsVisAttr = new G4VisAttributes(G4Colour::Magenta()); bumpsVisAttr->SetVisibility(true);

  // Epoxy betwen sensor and AlN cover
  G4double epoxySize_X = lgadSize_X;
  G4double epoxySize_Y = lgadSize_Y;
  G4double epoxySize_Z = 0.08*mm;
  G4Material* epoxyMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* epoxyVisAttr = new G4VisAttributes(G4Colour::Brown()); epoxyVisAttr->SetVisibility(true);

  // AlN sensor cover
  G4double coverplateSize_X = lgadSize_X;
  G4double coverplateSize_Y = lgadSize_Y; // Slightly wrong, should be slightly longer
  G4double coverplateSize_Z = 0.51*mm;
  G4Material* coverplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* coverplateVisAttr = new G4VisAttributes(G4Colour::Gray()); coverplateVisAttr->SetVisibility(true);

  G4double moduleSize_X = baseplateSize_X;
  G4double moduleSize_Y = baseplateSize_Y;
  G4double moduleSize_Z = basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+epoxySize_Z+lgadSize_Z+coverplateSize_Z;
  G4Material* moduleMat = G4Material::GetMaterial("Vacuum");
  G4VisAttributes* moduleVisAttr = new G4VisAttributes(G4Colour::Black()); moduleVisAttr->SetVisibility(false);

  // Build the module
  moduleBox = new G4Box(
    "ETLOneSensorModule",
    moduleSize_X/2., moduleSize_Y/2., moduleSize_Z/2.
  );
  moduleLogical = new G4LogicalVolume(
    moduleBox, 
    moduleMat,
    "ETLOneSensorModule"
  );
  moduleLogical->SetVisAttributes(moduleVisAttr);
  modulePV = new G4PVPlacement(
    0, relativePos,
    moduleLogical, "ETLOneSensorModule", motherLogical,
    false, 0, fCheckOverlaps
  );

  // Build the base film at the very bottom of the module
  G4Box* basefilmBox = new G4Box(
    "ETLOneSensorModuleBaseFilm",
    basefilmSize_X/2., basefilmSize_Y/2., basefilmSize_Z/2.
  );
  G4LogicalVolume* basefilmLogical = new G4LogicalVolume(
    basefilmBox,
    basefilmMat,
    "ETLOneSensorModuleBaseFilm"
  );
  basefilmLogical->SetVisAttributes(basefilmVisAttr);
  G4PVPlacement* basefilmPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, (basefilmSize_Z-moduleSize_Z)/2.),
    basefilmLogical, "ETLOneSensorModuleBaseFilm", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the base Al plate on top of the film
  G4Box* baseplateBox = new G4Box(
    "ETLOneSensorModuleBasePlate",
    baseplateSize_X/2., baseplateSize_Y/2., baseplateSize_Z/2.
  );
  G4LogicalVolume* baseplateLogical = new G4LogicalVolume(
    baseplateBox,
    baseplateMat,
    "ETLOneSensorModuleBasePlate"
  );
  baseplateLogical->SetVisAttributes(baseplateVisAttr);
  G4PVPlacement* baseplatePV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, basefilmSize_Z+(baseplateSize_Z-moduleSize_Z)/2.),
    baseplateLogical, "ETLOneSensorModuleBasePlate", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the laird film on top of the Al plate
  G4Box* lairdfilmBox = new G4Box(
    "ETLOneSensorModuleLairdFilm",
    lairdfilmSize_X/2., lairdfilmSize_Y/2., lairdfilmSize_Z/2.
  );
  G4LogicalVolume* lairdfilmLogical = new G4LogicalVolume(
    lairdfilmBox,
    lairdfilmMat,
    "ETLOneSensorModuleLairdFilm"
  );
  lairdfilmLogical->SetVisAttributes(lairdfilmVisAttr);
  G4PVPlacement* lairdfilmPV = new G4PVPlacement(
    0, G4ThreeVector(0, (lairdfilmSize_Y-baseplateSize_Y)/2., basefilmSize_Z+baseplateSize_Z+(lairdfilmSize_Z-moduleSize_Z)/2.),
    lairdfilmLogical, "ETLOneSensorModuleLairdFilm", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the ETROCs
  G4Box* etrocBox = new G4Box(
    "ETLOneSensorModuleETROC",
    etrocSize_X/2., etrocSize_Y/2., etrocSize_Z/2.
  );
  G4LogicalVolume* etrocLogical = new G4LogicalVolume(
    etrocBox,
    etrocMat,
    "ETLOneSensorModuleETROC"
  );
  etrocLogical->SetVisAttributes(etrocVisAttr);
  G4PVPlacement* etroc1PV = new G4PVPlacement(
    0, G4ThreeVector(-etrocSize_X/2., (etrocSize_Y-baseplateSize_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical, "ETLOneSensorModuleETROC1", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc2PV = new G4PVPlacement(
    0, G4ThreeVector(+etrocSize_X/2., (etrocSize_Y-baseplateSize_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical, "ETLOneSensorModuleETROC2", moduleLogical,
    false, 0, fCheckOverlaps
  );


  // Build the solder bumps
  G4Box* bumpsBox = new G4Box(
    "ETLOneSensorModuleSolderBumps",
    bumpsSize_X/2., bumpsSize_Y/2., bumpsSize_Z/2.
  );
  G4LogicalVolume* bumpsLogical = new G4LogicalVolume(
    bumpsBox,
    bumpsMat,
    "ETLOneSensorModuleSolderBumps"
  );
  bumpsLogical->SetVisAttributes(bumpsVisAttr);
  G4PVPlacement* bumpsPV = new G4PVPlacement(
    0, G4ThreeVector(0, (bumpsSize_Y-baseplateSize_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical, "ETLOneSensorModuleSolderBumps", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the LGAD
  G4Box* lgadBox = new G4Box(
    "ETLOneSensorModuleLGAD",
    lgadSize_X/2., lgadSize_Y/2., lgadSize_Z/2.
  );
  G4LogicalVolume* lgadLogical = new G4LogicalVolume(
    lgadBox,
    lgadMat,
    "ETLOneSensorModuleLGAD"
  );
  lgadLogical->SetVisAttributes(lgadVisAttr);
  G4PVPlacement* lgadPV = new G4PVPlacement(
    0, G4ThreeVector(0, (lgadSize_Y-baseplateSize_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
    lgadLogical, "ETLOneSensorModuleLGAD", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the epoxy cover
  G4Box* epoxyBox = new G4Box(
    "ETLOneSensorModuleEpoxyCover",
    epoxySize_X/2., epoxySize_Y/2., epoxySize_Z/2.
  );
  G4LogicalVolume* epoxyLogical = new G4LogicalVolume(
    epoxyBox,
    epoxyMat,
    "ETLOneSensorModuleEpoxyCover"
  );
  epoxyLogical->SetVisAttributes(epoxyVisAttr);
  G4PVPlacement* epoxyPV = new G4PVPlacement(
    0, G4ThreeVector(0, (epoxySize_Y-baseplateSize_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
    epoxyLogical, "ETLOneSensorModuleEpoxyCover", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the cover plate
  G4Box* coverplateBox = new G4Box(
    "ETLOneSensorModuleCoverPlate",
    coverplateSize_X/2., coverplateSize_Y/2., coverplateSize_Z/2.
  );
  G4LogicalVolume* coverplateLogical = new G4LogicalVolume(
    coverplateBox,
    coverplateMat,
    "ETLOneSensorModuleCoverPlate"
  );
  coverplateLogical->SetVisAttributes(coverplateVisAttr);
  G4PVPlacement* coverplatePV = new G4PVPlacement(
    0, G4ThreeVector(0, (coverplateSize_Y-baseplateSize_Y)/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical, "ETLOneSensorModuleCoverPlate", moduleLogical,
    false, 0, fCheckOverlaps
  );
}

void DetectorConstruction::BuildTwoSensorModule(
  G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos,
  G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV
){
  // See Fig. 3.56 in the MTD TDR for the module diagrams
  // Also Fig. 3.73 for the z positions.
  // Note that the base aluminum plate on which the service is loaded, or the epoxy underneath are not included since they are only on one side.
  // AlN
  G4double baseplateSize_X = 43.1*mm;
  G4double baseplateSize_Y = 56.5*mm;
  G4double baseplateSize_Z = 0.79*mm;
  G4Material* baseplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* baseplateVisAttr = new G4VisAttributes(G4Colour::Gray()); baseplateVisAttr->SetVisibility(true);

  // Thermal pad
  G4double basefilmSize_X = baseplateSize_X;
  G4double basefilmSize_Y = baseplateSize_Y;
  G4double basefilmSize_Z = 0.25*mm;
  G4Material* basefilmMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* basefilmVisAttr = new G4VisAttributes(G4Colour::Brown()); basefilmVisAttr->SetVisibility(true);

  // ETROC
  G4double etrocSize_X = 20.8*mm;
  G4double etrocSize_Y = 22.3*mm;
  G4double etrocSize_Z = 0.25*mm;
  G4Material* etrocMat = G4Material::GetMaterial("Si");
  G4VisAttributes* etrocVisAttr = new G4VisAttributes(G4Colour::Green()); etrocVisAttr->SetVisibility(true);

  // Laird film
  G4double lairdfilmSize_X = etrocSize_X*2.;
  G4double lairdfilmSize_Y = etrocSize_Y*2.;
  G4double lairdfilmSize_Z = 0.08*mm;
  G4Material* lairdfilmMat = G4Material::GetMaterial("Laird");
  G4VisAttributes* lairdfilmVisAttr = new G4VisAttributes(G4Colour::Brown()); lairdfilmVisAttr->SetVisibility(true);

  // LGAD sensor
  G4double lgadSize_X = 42.0*mm;
  G4double lgadSize_Y = 21.2*mm;
  G4double lgadSize_Z = 0.3*mm;
  G4Material* lgadMat = G4Material::GetMaterial("Si");
  G4VisAttributes* lgadVisAttr = new G4VisAttributes(G4Colour::Yellow()); lgadVisAttr->SetVisibility(true);

  // Solder bumps
  G4double bumpsSize_X = etrocSize_X*2.;
  G4double bumpsSize_Y = lgadSize_Y;
  G4double bumpsSize_Z = 0.03*mm;
  G4Material* bumpsMat = G4Material::GetMaterial("SensorBump");
  G4VisAttributes* bumpsVisAttr = new G4VisAttributes(G4Colour::Magenta()); bumpsVisAttr->SetVisibility(true);

  // Epoxy betwen sensor and AlN cover
  G4double epoxySize_X = lgadSize_X;
  G4double epoxySize_Y = lgadSize_Y;
  G4double epoxySize_Z = 0.08*mm;
  G4Material* epoxyMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* epoxyVisAttr = new G4VisAttributes(G4Colour::Brown()); epoxyVisAttr->SetVisibility(true);

  // AlN sensor cover
  G4double coverplateSize_X = lgadSize_X;
  G4double coverplateSize_Y = lgadSize_Y; // Slightly wrong, should be slightly longer
  G4double coverplateSize_Z = 0.51*mm;
  G4Material* coverplateMat = G4Material::GetMaterial("AlN");
  G4VisAttributes* coverplateVisAttr = new G4VisAttributes(G4Colour::Gray()); coverplateVisAttr->SetVisibility(true);

  G4double moduleSize_X = baseplateSize_X;
  G4double moduleSize_Y = baseplateSize_Y;
  G4double moduleSize_Z = basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+epoxySize_Z+lgadSize_Z+coverplateSize_Z;
  G4Material* moduleMat = G4Material::GetMaterial("Vacuum");
  G4VisAttributes* moduleVisAttr = new G4VisAttributes(G4Colour::Black()); moduleVisAttr->SetVisibility(false);

  // Build the module
  moduleBox = new G4Box(
    "ETLTwoSensorModule",
    moduleSize_X/2., moduleSize_Y/2., moduleSize_Z/2.
  );
  moduleLogical = new G4LogicalVolume(
    moduleBox,
    moduleMat,
    "ETLTwoSensorModule"
  );
  moduleLogical->SetVisAttributes(moduleVisAttr);
  modulePV = new G4PVPlacement(
    0, relativePos,
    moduleLogical, "ETLTwoSensorModule", motherLogical,
    false, 0, fCheckOverlaps
  );

  // Build the base film at the very bottom of the module
  G4Box* basefilmBox = new G4Box(
    "ETLTwoSensorModuleBaseFilm",
    basefilmSize_X/2., basefilmSize_Y/2., basefilmSize_Z/2.
  );
  G4LogicalVolume* basefilmLogical = new G4LogicalVolume(
    basefilmBox,
    basefilmMat,
    "ETLTwoSensorModuleBaseFilm"
  );
  basefilmLogical->SetVisAttributes(basefilmVisAttr);
  G4PVPlacement* basefilmPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, (basefilmSize_Z-moduleSize_Z)/2.),
    basefilmLogical, "ETLTwoSensorModuleBaseFilm", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the base Al plate on top of the film
  G4Box* baseplateBox = new G4Box(
    "ETLTwoSensorModuleBasePlate",
    baseplateSize_X/2., baseplateSize_Y/2., baseplateSize_Z/2.
  );
  G4LogicalVolume* baseplateLogical = new G4LogicalVolume(
    baseplateBox,
    baseplateMat,
    "ETLTwoSensorModuleBasePlate"
  );
  baseplateLogical->SetVisAttributes(baseplateVisAttr);
  G4PVPlacement* baseplatePV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, basefilmSize_Z+(baseplateSize_Z-moduleSize_Z)/2.),
    baseplateLogical, "ETLTwoSensorModuleBasePlate", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the laird films on top of the Al plate
  G4Box* lairdfilmBox = new G4Box(
    "ETLTwoSensorModuleLairdFilm",
    lairdfilmSize_X/2., lairdfilmSize_Y/2., lairdfilmSize_Z/2.
  );
  G4LogicalVolume* lairdfilmLogical = new G4LogicalVolume(
    lairdfilmBox,
    lairdfilmMat,
    "ETLTwoSensorModuleLairdFilm"
  );
  lairdfilmLogical->SetVisAttributes(lairdfilmVisAttr);
  G4PVPlacement* lairdfilmPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, basefilmSize_Z+baseplateSize_Z+(lairdfilmSize_Z-moduleSize_Z)/2.),
    lairdfilmLogical, "ETLTwoSensorModuleLairdFilm", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the ETROCs
  G4Box* etrocBox = new G4Box(
    "ETLTwoSensorModuleETROC",
    etrocSize_X/2., etrocSize_Y/2., etrocSize_Z/2.
  );
  G4LogicalVolume* etrocLogical = new G4LogicalVolume(
    etrocBox,
    etrocMat,
    "ETLTwoSensorModuleETROC"
  );
  etrocLogical->SetVisAttributes(etrocVisAttr);
  G4PVPlacement* etroc1PV = new G4PVPlacement(
    0, G4ThreeVector(-etrocSize_X/2., -etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical, "ETLTwoSensorModuleETROC1", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc2PV = new G4PVPlacement(
    0, G4ThreeVector(+etrocSize_X/2., -etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical, "ETLTwoSensorModuleETROC2", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc3PV = new G4PVPlacement(
    0, G4ThreeVector(-etrocSize_X/2., etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical, "ETLTwoSensorModuleETROC3", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* etroc4PV = new G4PVPlacement(
    0, G4ThreeVector(+etrocSize_X/2., etrocSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+(etrocSize_Z-moduleSize_Z)/2.),
    etrocLogical, "ETLTwoSensorModuleETROC4", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the solder bumps
  G4Box* bumpsBox = new G4Box(
    "ETLTwoSensorModuleSolderBumps",
    bumpsSize_X/2., bumpsSize_Y/2., bumpsSize_Z/2.
  );
  G4LogicalVolume* bumpsLogical = new G4LogicalVolume(
    bumpsBox,
    bumpsMat,
    "ETLTwoSensorModuleSolderBumps"
  );
  bumpsLogical->SetVisAttributes(bumpsVisAttr);
  G4PVPlacement* bumps1PV = new G4PVPlacement(
    0, G4ThreeVector(0, -bumpsSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical, "ETLTwoSensorModuleSolderBumps1", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* bumps2PV = new G4PVPlacement(
    0, G4ThreeVector(0, +bumpsSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+(bumpsSize_Z-moduleSize_Z)/2.),
    bumpsLogical, "ETLTwoSensorModuleSolderBumps2", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the LGAD
  G4Box* lgadBox = new G4Box(
    "ETLTwoSensorModuleLGAD",
    lgadSize_X/2., lgadSize_Y/2., lgadSize_Z/2.
  );
  G4LogicalVolume* lgadLogical = new G4LogicalVolume(
    lgadBox,
    lgadMat,
    "ETLTwoSensorModuleLGAD"
  );
  lgadLogical->SetVisAttributes(lgadVisAttr);
  G4PVPlacement* lgad1PV = new G4PVPlacement(
    0, G4ThreeVector(0, -lgadSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
    lgadLogical, "ETLTwoSensorModuleLGAD1", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* lgad2PV = new G4PVPlacement(
    0, G4ThreeVector(0, +lgadSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+(lgadSize_Z-moduleSize_Z)/2.),
    lgadLogical, "ETLTwoSensorModuleLGAD2", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the epoxy cover
  G4Box* epoxyBox = new G4Box(
    "ETLTwoSensorModuleEpoxyCover",
    epoxySize_X/2., epoxySize_Y/2., epoxySize_Z/2.
  );
  G4LogicalVolume* epoxyLogical = new G4LogicalVolume(
    epoxyBox,
    epoxyMat,
    "ETLTwoSensorModuleEpoxyCover"
  );
  epoxyLogical->SetVisAttributes(epoxyVisAttr);
  G4PVPlacement* epoxy1PV = new G4PVPlacement(
    0, G4ThreeVector(0, -epoxySize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
    epoxyLogical, "ETLTwoSensorModuleEpoxyCover1", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* epoxy2PV = new G4PVPlacement(
    0, G4ThreeVector(0, +epoxySize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+(epoxySize_Z-moduleSize_Z)/2.),
    epoxyLogical, "ETLTwoSensorModuleEpoxyCover2", moduleLogical,
    false, 0, fCheckOverlaps
  );

  // Build the cover plate
  G4Box* coverplateBox = new G4Box(
    "ETLTwoSensorModuleCoverPlate",
    coverplateSize_X/2., coverplateSize_Y/2., coverplateSize_Z/2.
  );
  G4LogicalVolume* coverplateLogical = new G4LogicalVolume(
    coverplateBox,
    coverplateMat,
    "ETLTwoSensorModuleCoverPlate"
  );
  coverplateLogical->SetVisAttributes(coverplateVisAttr);
  G4PVPlacement* coverplate1PV = new G4PVPlacement(
    0, G4ThreeVector(0, -coverplateSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical, "ETLTwoSensorModuleCoverPlate1", moduleLogical,
    false, 0, fCheckOverlaps
  );
  G4PVPlacement* coverplate2PV = new G4PVPlacement(
    0, G4ThreeVector(0, +coverplateSize_Y/2., basefilmSize_Z+baseplateSize_Z+lairdfilmSize_Z+etrocSize_Z+bumpsSize_Z+lgadSize_Z+epoxySize_Z+(coverplateSize_Z-moduleSize_Z)/2.),
    coverplateLogical, "ETLTwoSensorModuleCoverPlate2", moduleLogical,
    false, 0, fCheckOverlaps
  );
}


void DetectorConstruction::BuildSensorServiceHybrid(
  int const& nSensorsPerSide, // 6 or 3
  G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos,
  G4Box*& serviceBox, G4LogicalVolume*& serviceLogical, G4PVPlacement*& servicePV
){
  // See Figs. 3.61 and 3.62 in the MTD TDR for the diagrams
  G4double servicehybridSize_X = 56.5*mm; // FIXME: Not sure what the width of the service hybrid is
  G4double servicehybridSize_Y = 43.1*((G4double) nSensorsPerSide)*mm;

  // Thermal pad
  G4double thermalpadSize_X = servicehybridSize_X;
  G4double thermalpadSize_Y = servicehybridSize_Y;
  G4double thermalpadSize_Z = 0.25*mm;
  G4Material* thermalpadMat = G4Material::GetMaterial("Epoxy");
  G4VisAttributes* thermalpadVisAttr = new G4VisAttributes(G4Colour::Brown()); thermalpadVisAttr->SetVisibility(true);

  // Readout board
  G4double readoutboardSize_X = servicehybridSize_X;
  G4double readoutboardSize_Y = servicehybridSize_Y;
  G4double readoutboardSize_Z = 1.*mm;
  G4Material* readoutboardMat = G4Material::GetMaterial("PCB");
  G4VisAttributes* readoutboardVisAttr = new G4VisAttributes(G4Colour::Brown()); readoutboardVisAttr->SetVisibility(true);

  // Connectors/flex/stiffener
  G4double connectorSize_X = servicehybridSize_X;
  G4double connectorSize_Y = servicehybridSize_Y;
  G4double connectorSize_Z = 1.5*mm;
  G4Material* connectorMat = G4Material::GetMaterial("ServiceConnector");
  G4VisAttributes* connectorVisAttr = new G4VisAttributes(G4Colour::Brown()); connectorVisAttr->SetVisibility(true);

  // Power board
  G4double powerboardSize_X = servicehybridSize_X;
  G4double powerboardSize_Y = servicehybridSize_Y;
  G4double powerboardSize_Z = 3.1*mm;
  G4Material* powerboardMat = G4Material::GetMaterial("PCB");
  G4VisAttributes* powerboardVisAttr = new G4VisAttributes(G4Colour::Brown()); powerboardVisAttr->SetVisibility(true);

  G4double servicehybridSize_Z = (thermalpadSize_Z+readoutboardSize_Z+connectorSize_Z+powerboardSize_Z);
  G4Material* servicehybridMat = G4Material::GetMaterial("Vacuum");
  G4VisAttributes* servicehybridVisAttr = new G4VisAttributes(G4Colour::Gray()); servicehybridVisAttr->SetVisibility(true);

  // Build the module
  serviceBox = new G4Box(
    "ETL12SensorServiceHybrid",
    servicehybridSize_X/2., servicehybridSize_Y/2., servicehybridSize_Z/2.
  );
  serviceLogical = new G4LogicalVolume(
    serviceBox,
    servicehybridMat,
    "ETL12SensorServiceHybrid"
  );
  serviceLogical->SetVisAttributes(servicehybridVisAttr);
  servicePV = new G4PVPlacement(
    0, relativePos,
    serviceLogical, "ETL12SensorServiceHybrid", motherLogical,
    false, 0, fCheckOverlaps
  );


  // Thermal pad
  G4Box* thermalpadBox = new G4Box(
    "ETL12SensorServiceHybrid_ThermalPad",
    thermalpadSize_X/2., thermalpadSize_Y/2., thermalpadSize_Z/2.
  );
  G4LogicalVolume* thermalpadLogical = new G4LogicalVolume(
    thermalpadBox,
    thermalpadMat,
    "ETL12SensorServiceHybrid_ThermalPad"
  );
  thermalpadLogical->SetVisAttributes(thermalpadVisAttr);
  G4PVPlacement* thermalpadPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, (thermalpadSize_Z-serviceSize_Y)/2.),
    thermalpadLogical, "ETL12SensorServiceHybrid_ThermalPad", serviceLogical,
    false, 0, fCheckOverlaps
  );

  // Readout board
  G4Box* readoutboardBox = new G4Box(
    "ETL12SensorServiceHybrid_ReadoutBoard",
    readoutboardSize_X/2., readoutboardSize_Y/2., readoutboardSize_Z/2.
  );
  G4LogicalVolume* readoutboardLogical = new G4LogicalVolume(
    readoutboardBox,
    readoutboardMat,
    "ETL12SensorServiceHybrid_ReadoutBoard"
  );
  readoutboardLogical->SetVisAttributes(readoutboardVisAttr);
  G4PVPlacement* readoutboardPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, thermalpadSize_Z+(readoutboardSize_Z-serviceSize_Y)/2.),
    readoutboardLogical, "ETL12SensorServiceHybrid_ReadoutBoard", serviceLogical,
    false, 0, fCheckOverlaps
  );

  // Connectors
  G4Box* connectorBox = new G4Box(
    "ETL12SensorServiceHybrid_Connectors",
    connectorSize_X/2., connectorSize_Y/2., connectorSize_Z/2.
  );
  G4LogicalVolume* connectorLogical = new G4LogicalVolume(
    connectorBox,
    connectorMat,
    "ETL12SensorServiceHybrid_Connectors"
  );
  connectorLogical->SetVisAttributes(connectorVisAttr);
  G4PVPlacement* connectorPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, thermalpadSize_Z+readoutboardSize_Z+(connectorSize_Z-serviceSize_Y)/2.),
    connectorLogical, "ETL12SensorServiceHybrid_Connectors", serviceLogical,
    false, 0, fCheckOverlaps
  );

  // Power board
  G4Box* powerboardBox = new G4Box(
    "ETL12SensorServiceHybrid_PowerBoard",
    powerboardSize_X/2., powerboardSize_Y/2., powerboardSize_Z/2.
  );
  G4LogicalVolume* powerboardLogical = new G4LogicalVolume(
    powerboardBox,
    powerboardMat,
    "ETL12SensorServiceHybrid_PowerBoard"
  );
  powerboardLogical->SetVisAttributes(powerboardVisAttr);
  G4PVPlacement* powerboardPV = new G4PVPlacement(
    0, G4ThreeVector(0, 0, thermalpadSize_Z+readoutboardSize_Z+connectorSize_Z+(powerboardSize_Z-serviceSize_Y)/2.),
    powerboardLogical, "ETL12SensorServiceHybrid_PowerBoard", serviceLogical,
    false, 0, fCheckOverlaps
  );
}



G4VPhysicalVolume* DetectorConstruction::DefineVolumes(){
  // Setup materials
  G4Material* world_mat = G4Material::GetMaterial("Vacuum");
  G4Material* env_mat = G4Material::GetMaterial("G4_CARBON_DIOXIDE");
  if (!world_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve world_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }
  if (!env_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve env_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }

  // Size parameters
  G4double env_sizeXY = 5*cm, env_sizeZ = 0.3*cm;
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;

  // World
  G4Box* solidWorld = new G4Box(
    "World",
    0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ
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
  G4Box* solidEnv = new G4Box(
    "Envelope",                                   //its name
    0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ //its size
  );
  G4LogicalVolume* logicEnv = new G4LogicalVolume(
    solidEnv,            //its solid
    env_mat,             //its material
    "Envelope"           //its name
  );
  G4VisAttributes* envVisAttr = new G4VisAttributes(G4Colour::Black()); envVisAttr->SetVisibility(false);
  logicEnv->SetVisAttributes(envVisAttr);
  new G4PVPlacement(
    0,                       //no rotation
    G4ThreeVector(),         //at (0,0,0)
    logicEnv,                //its logical volume
    "Envelope",              //its name
    logicWorld,              //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    fCheckOverlaps            //overlaps checking
  );

  G4Box* servicehybrid6Box;
  G4LogicalVolume* servicehybrid6Logical;
  G4PVPlacement* servicehybrid6PV;
  BuildSensorServiceHybrid(
    6,
    logicEnv, G4ThreeVector(0, 0, 0),
    servicehybrid6Box, servicehybrid6Logical, servicehybrid6PV
  );

  /*
  G4Box* servicehybrid3Box;
  G4LogicalVolume* servicehybrid3Logical;
  G4PVPlacement* servicehybrid3PV;
  BuildSensorServiceHybrid(
    3,
    logicEnv, G4ThreeVector(0, 0, 0),
    servicehybrid3Box, servicehybrid3Logical, servicehybrid3PV
  );
  */

  /*
  G4Box* moduleBox;
  G4LogicalVolume* moduleLogical;
  G4PVPlacement* modulePV;
  BuildOneSensorModule(
    logicEnv, G4ThreeVector(0,0,0),
    moduleBox, moduleLogical, modulePV
  );
  */

  // Set Shape2 as scoring volume
  fScoringVolume = moduleLogical;

  //always return the physical World
  return physWorld;
}
