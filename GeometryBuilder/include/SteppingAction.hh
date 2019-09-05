#ifndef SteppingAction_h
#define SteppingAction_h

#include <vector>
#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class G4LogicalVolume;

/// Stepping action class
class SteppingAction : public G4UserSteppingAction{
public:
  SteppingAction(EventAction* eventAction);
  virtual ~SteppingAction();

  // method from the base class
  virtual void UserSteppingAction(const G4Step*);

private:
  EventAction* fEventAction;
  std::vector<G4LogicalVolume*> fScoringVolumes;
};


#endif
