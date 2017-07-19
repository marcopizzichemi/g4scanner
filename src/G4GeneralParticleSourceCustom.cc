#include "G4GeneralParticleSourceCustom.hh"

void G4GeneralParticleSourceCustom::GeneratePrimaryVertex(G4Event* evt)
{
     GetCurrentSource()->GeneratePrimaryVertex(evt);
}
