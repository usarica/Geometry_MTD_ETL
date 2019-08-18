#ifndef DetectorConstruction_h
#define DetectorConstruction_h

#include <vector>
#include <utility>
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"


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
  void BuildOneSensorModule(bool rightFlank, G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos, G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV);
  void BuildTwoSensorModule(G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos, G4Box*& moduleBox, G4LogicalVolume*& moduleLogical, G4PVPlacement*& modulePV);
  void BuildSensorServiceHybrid(int const& nSensorsPerSide, G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos, G4Box*& serviceBox, G4LogicalVolume*& serviceLogical, G4PVPlacement*& servicePV);
  void BuildWedgeComponents(G4LogicalVolume* motherLogical, std::vector<std::pair<G4double, G4double>> const& coolingpipes_xpos_ymin);

public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  G4LogicalVolume* GetScoringVolume() const{ return fScoringVolume; }

};

#endif
