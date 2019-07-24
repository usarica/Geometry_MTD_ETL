#ifndef DetectorConstruction_h
#define DetectorConstruction_h

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;


// Detector construction class to define materials and geometry.
class DetectorConstruction : public G4VUserDetectorConstruction{
protected:
  G4bool fCheckOverlaps; // Option to switch on/off checking of volumes overlaps
  G4LogicalVolume* fScoringVolume;

  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();
  void BuildOneSensorModule(G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos, G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV);
  void BuildTwoSensorModule(G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos, G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV);
  void BuildSensorServiceHybrid(int const& nSensorsPerSide, G4LogicalVolume* motherLogical, G4ThreeVector const& relativePos, G4Box*& serviceBox, G4LogicalVolume*& serviceLogical, G4PVPlacement*& servicePV);

public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  G4LogicalVolume* GetScoringVolume() const{ return fScoringVolume; }

};

#endif
