#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"


class RunAction;


// Event action class
class EventAction : public G4UserEventAction{
private:
  RunAction* fRunAction;
  G4double fEdep;

public:
  EventAction(RunAction* runAction);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  void AddEdep(G4double edep){ fEdep += edep; }

};


#endif
