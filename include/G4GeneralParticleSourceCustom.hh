#include "G4GeneralParticleSource.hh"

class G4GeneralParticleSourceCustom : public G4GeneralParticleSource
{
public:
  void GeneratePrimaryVertex(G4Event* evt);
};
