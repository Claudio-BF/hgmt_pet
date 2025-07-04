
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

#ifndef AnnihilScorer_hh
#define AnnihilScorer_hh

#include "G4Types.hh"
#include "TsVNtupleScorer.hh"

#include <cmath>
#include <list>
using namespace std;
struct vec3d {
  double x, y, z;
};
struct event {
  double tof;
  double energy;
  vec3d position;
  vec3d direction;
  uint primary;
};
struct prim_annihilation {
  double time_of_annihilation;
  vec3d origin;
  vec3d center;
  vector<event> events;
};
class AnnihilScorer : public TsVNtupleScorer {
public:
  AnnihilScorer(TsParameterManager *pM, TsMaterialManager *mM,
                TsGeometryManager *gM, TsScoringManager *scM,
                TsExtensionManager *eM, G4String scorerName, G4String quantity,
                G4String outFileName, G4bool isSubScorer);

  virtual ~AnnihilScorer();

  G4bool ProcessHits(G4Step *, G4TouchableHistory *);
  // called after the last hit of a given track
  void UserHookForEndOfTrack(const G4Track *);
  // called after the last hit of all tracks resulting from a given particle
  // incident on the scoring component
  void UserHookForEndOfIncidentParticle();
  // called at the end of the event
  void UserHookForEndOfEvent();
  // called at the end of the run
  void UserHookForEndOfRun();

private:
  prim_annihilation new_annihilation;
  int trackid1 = -1, trackid2 = -1;
};
#endif
