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

  nist->FindOrBuildMaterial("G4_AIR");
  nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  G4double env_temp = NTP_Temperature; // 293.15*kelvin
  G4double env_pressure = CLHEP::STP_Pressure; // in units of pascal

  // Vacuum
  new G4Material(
    "Galactic", z=1., a=1.01*g/mole, density=universe_mean_density,
    kStateGas, 2.73*kelvin, 3.e-18*pascal
  );

  G4cout << "DetectorConstruction::DefineMaterials: Environment materials built." << G4endl;

  // MIC6 Aluminum
  // Used as ETL disk material
  nist->BuildMaterialWithNewDensity("MIC6_Al", "G4_Al", 2.796*g/cm3, env_temp, env_pressure); // 0.101*lb/in3

  // When calling FindOrBuildElement, do not use 'G4_' at the beginning of the element name
  G4Element* Aluminum = nist->FindOrBuildElement("Al");
  G4Element* Nitrogen = nist->FindOrBuildElement("N");
  G4cout << "DetectorConstruction::DefineMaterials: Elements built." << G4endl;

  if (!Aluminum) G4cerr << "DetectorConstruction::DefineMaterials: Element Al not found!" << G4endl;
  if (!Nitrogen) G4cerr << "DetectorConstruction::DefineMaterials: Element N not found!" << G4endl;

  // Aluminum Nitride (AlN), used in sensors and elsewhere
  G4Material* AlN = new G4Material("AlN", 3.260*g/cm3, 2, kStateSolid, env_temp, env_pressure);
  AlN->AddElement(Aluminum, (G4int) 1);
  AlN->AddElement(Nitrogen, (G4int) 1);

  G4cout << "DetectorConstruction::DefineMaterials: AlN built!" << G4endl;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes(){
  // Setup materials
  G4Material* world_mat = G4Material::GetMaterial("Galactic");
  G4Material* env_mat = G4Material::GetMaterial("G4_CARBON_DIOXIDE");
  G4Material* shape1_mat = G4Material::GetMaterial("MIC6_Al");
  G4Material* shape2_mat = G4Material::GetMaterial("AlN");
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
  if (!shape1_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve shape1_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }
  if (!shape2_mat){
    G4ExceptionDescription msg;
    msg << "Cannot retrieve shape2_mat.";
    G4Exception("DetectorConstruction::DefineVolumes",
                "MyCode0001", FatalException, msg);
  }

  // Size parameters
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
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

  // Shape 1
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);

  // Conical section shape       
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Cons* solidShape1 = new G4Cons(
    "Shape1",
    shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
    shape1_phimin, shape1_phimax
  );
  G4LogicalVolume* logicShape1 = new G4LogicalVolume(
    solidShape1,         //its solid
    shape1_mat,          //its material
    "Shape1"             //its name
  );
  new G4PVPlacement(
    0,                       //no rotation
    pos1,                    //at position
    logicShape1,             //its logical volume
    "Shape1",                //its name
    logicEnv,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    fCheckOverlaps            //overlaps checking
  );

  // Shape 2
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape       
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;
  G4Trd* solidShape2 = new G4Trd(
    "Shape2",                                       //its name
    0.5*shape2_dxa, 0.5*shape2_dxb,
    0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz   //its size
  );
  G4LogicalVolume* logicShape2 = new G4LogicalVolume(
    solidShape2,         //its solid
    shape2_mat,          //its material
    "Shape2"             //its name
  );
  new G4PVPlacement(
    0,                       //no rotation
    pos2,                    //at position
    logicShape2,             //its logical volume
    "Shape2",                //its name
    logicEnv,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    fCheckOverlaps            //overlaps checking
  );

  // Set Shape2 as scoring volume
  fScoringVolume = logicShape2;

  //always return the physical World
  return physWorld;
}
