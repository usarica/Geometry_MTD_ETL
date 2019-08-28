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
private:
  constexpr static bool putServiceHybrids = true;
  constexpr static bool putModules = true;
  constexpr static bool doWedgeFarFace = true;
  constexpr static bool doWedgeCloseFace = true;
  constexpr static bool doFullDisk = true;
  constexpr static bool doFarDisk = true;
  constexpr static bool doCloseDisk = true;
  constexpr static bool putWedgeComponents = true;
  constexpr static bool addMIC6Al = true;
  constexpr static bool addEpoxy = true;
  constexpr static bool addCoolingAl = true;
  constexpr static bool drillCoolingPipeCavities = true;
  constexpr static bool addCoolingPipes = true && DetectorConstruction::drillCoolingPipeCavities;
  constexpr static bool addFrontFaceSupportBars = true && DetectorConstruction::doWedgeFarFace;
  constexpr static bool doEndCapSupport = true;
  constexpr static bool addSupportTube = true && DetectorConstruction::doEndCapSupport;
  constexpr static bool addMountingBracket = true && DetectorConstruction::doEndCapSupport;
  constexpr static bool doFrontEndcap = true;
  constexpr static bool doBackEndcap = true;
  constexpr static bool endcapsAreReflected = true;

protected:
  G4bool fCheckOverlaps; // Option to switch on/off checking of volumes overlaps
  G4LogicalVolume* fScoringVolume;

  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();
  void BuildOneSensorModule(G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos);
  void BuildTwoSensorModule(G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos);
  void BuildSensorServiceHybrid(int const& nSensorsPerSide, G4LogicalVolume* motherLogical, G4RotationMatrix* rotation, G4ThreeVector const& relativePos);
  void BuildWedgeComponents(G4LogicalVolume* passiveLogical, G4LogicalVolume* activeFarLogical, std::vector<std::pair<G4double, G4double>> const& coolingpipes_xpos_ymin);
  void BuildEndCapSupportComponents(G4LogicalVolume* motherLogical);

public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  G4LogicalVolume* GetScoringVolume() const{ return fScoringVolume; }

};

#endif
