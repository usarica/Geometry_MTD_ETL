#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"


SteppingAction::SteppingAction(EventAction* eventAction) :
  G4UserSteppingAction(),
  fEventAction(eventAction)
{}
SteppingAction::~SteppingAction(){}

void SteppingAction::UserSteppingAction(const G4Step* step){
  if (fScoringVolumes.empty()){
    const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolumes = detectorConstruction->GetScoringVolumes();
  }

  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // check if we are in a scoring volume
  G4bool inVol = false;
  for (auto const& v:fScoringVolumes) inVol |= (volume == v);
  if (!inVol){
    G4cout << "Scoring volume " << volume->GetName() << " could not be found!" << G4endl;
    return;
  }
  else G4cout << "Scoring volume " << volume->GetName() << " is found!" << G4endl;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
}
