
// Scorer for AnnihilScorer
//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#include "AnnihilScorer.hh"
#include "G4ThreeVector.hh"
#include "G4VProcess.hh"
#include "TsTrackInformation.hh"
#include "TsVNtupleScorer.hh"
#include <cstdio>
#include <unistd.h>
using namespace std;
FILE *outFile;
prim_annihilation new_annihilation;
int trackid1 = -1, trackid2 = -1;
void write_vec3d(const vec3d &vec, FILE *file) {
  fwrite(&vec.x, sizeof(double), 1, file);
  fwrite(&vec.y, sizeof(double), 1, file);
  fwrite(&vec.z, sizeof(double), 1, file);
}
void write_event(const event &new_event, FILE *file) {
  fwrite(&new_event.tof, sizeof(double), 1, file);
  fwrite(&new_event.energy, sizeof(double), 1, file);
  write_vec3d(new_event.position, file);
  write_vec3d(new_event.direction, file);
  fwrite(&new_event.primary, sizeof(uint), 1, file);
}
void write_annihilation(const prim_annihilation &annihilation, FILE *file) {
  // Write fixed-size components
  fwrite(&annihilation.time_of_annihilation, sizeof(double), 1, file);
  write_vec3d(new_annihilation.origin, file);
  write_vec3d(new_annihilation.center, file);
  uint num_events = (uint)new_annihilation.events.size();
  fwrite(&num_events, sizeof(uint), 1, file);
  for (event new_event : annihilation.events)
    write_event(new_event, file);
  // Write each element individually
}
AnnihilScorer::AnnihilScorer(TsParameterManager *pM, TsMaterialManager *mM,
                             TsGeometryManager *gM, TsScoringManager *scM,
                             TsExtensionManager *eM, G4String scorerName,
                             G4String quantity, G4String outFileName,
                             G4bool isSubScorer)
    : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName,
                      isSubScorer) {
  outFile = fopen(outFileName, "wb");
  SuppressStandardOutputHandling();
}

AnnihilScorer::~AnnihilScorer() { fclose(outFile); }

G4bool AnnihilScorer::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  G4Track *track = aStep->GetTrack();
  G4TrackStatus trackStatus = track->GetTrackStatus();
  G4int particleType = track->GetDefinition()->GetPDGEncoding();
  G4int parentID = track->GetParentID();

  G4int stepNumber = track->GetCurrentStepNumber();
  // first gamma scattering
  if (particleType == 22 && stepNumber == 1 &&
      track->GetCreatorProcess()->GetProcessName() == "annihil") {
    if (trackid1 == -1)
      trackid1 = track->GetTrackID();
    else
      trackid2 = track->GetTrackID();
  }
  // electron first step in particle
  else if (stepNumber == 1 && particleType == 11) {
    event new_event;
    new_event.tof = track->GetGlobalTime();
    new_event.energy = track->GetVertexKineticEnergy() * 1000;
    G4ThreeVector threevec;
    threevec = track->GetVertexPosition();
    new_event.position.x = threevec.x() * 0.1;
    new_event.position.y = threevec.y() * 0.1;
    new_event.position.z = threevec.z() * 0.1;
    threevec = track->GetVertexMomentumDirection();
    new_event.direction.x = threevec.x();
    new_event.direction.y = threevec.y();
    new_event.direction.z = threevec.z();
    if (track->GetParentID() == trackid1)
      new_event.primary = 1;
    else if (track->GetParentID() == trackid2)
      new_event.primary = 2;
    else
      new_event.primary = 0;
    new_annihilation.events.push_back(new_event);
  }
  // final positron annihilation
  else if (particleType == -11 && parentID == 0 &&
           (trackStatus == fStopAndKill ||
            (trackStatus == fStopButAlive &&
             aStep->GetPostStepPoint()->GetKineticEnergy() == 0))) {
    new_annihilation.time_of_annihilation =
        aStep->GetPostStepPoint()->GetGlobalTime();
    G4ThreeVector threevec = aStep->GetPostStepPoint()->GetPosition();
    new_annihilation.center.x = threevec.x() * 0.1;
    new_annihilation.center.y = threevec.y() * 0.1;
    new_annihilation.center.z = threevec.z() * 0.1;
  }

  // initial positron
  else if (particleType == -11 && parentID == 0 && stepNumber == 1) {
    G4ThreeVector threevec = track->GetVertexPosition();
    new_annihilation.origin.x = threevec.x() * 0.1;
    new_annihilation.origin.y = threevec.y() * 0.1;
    new_annihilation.origin.z = threevec.z() * 0.1;
  }
  return true;
}
void AnnihilScorer::UserHookForEndOfTrack(const G4Track *) {}
void AnnihilScorer::UserHookForEndOfIncidentParticle() {}
void AnnihilScorer::UserHookForEndOfEvent() {
  sort(new_annihilation.events.begin(), new_annihilation.events.end(),
       [](const event &a, const event &b) { return a.tof < b.tof; });
  write_annihilation(new_annihilation, outFile);
  new_annihilation = prim_annihilation{};
  trackid1 = -1;
  trackid2 = -1;
}
void AnnihilScorer::UserHookForEndOfRun() {}
