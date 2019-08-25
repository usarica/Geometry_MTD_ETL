#ifndef DetectorConstruction_h
#define DetectorConstruction_h

#include <vector>
#include <utility>
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PVPlacement;


namespace DetectorConstructionHelpers{
  struct BasicDetectorAttributes{
    G4VSolid* solid;
    G4LogicalVolume* logical;
    G4Material* material;
    G4VisAttributes* visualization;

    BasicDetectorAttributes();
    BasicDetectorAttributes(BasicDetectorAttributes const& other);
    BasicDetectorAttributes(G4VSolid* solid_, G4LogicalVolume* logical_, G4Material* material_, G4VisAttributes* visualization_);
  };
}

// Detector construction class to define materials and geometry.
class DetectorConstruction : public G4VUserDetectorConstruction{
protected:
  G4bool fCheckOverlaps; // Option to switch on/off checking of volumes overlaps
  G4LogicalVolume* fScoringVolume;

  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();
  void BuildOneSensorModule(G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos);
  void BuildTwoSensorModule(G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos);
  void BuildSensorServiceHybrid(int const& nSensorsPerSide, G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos);
  void BuildWedgeComponents(G4LogicalVolume* passiveLogical, G4LogicalVolume* activeFarLogical, std::vector<std::pair<G4double, G4double>> const& coolingpipes_xpos_ymin);

public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  G4LogicalVolume* GetScoringVolume() const{ return fScoringVolume; }

};

#endif
