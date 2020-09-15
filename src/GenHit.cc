#include "GenHit.hh"

G4ThreadLocal G4Allocator<GenHit>* GenHitAllocator;

GenHit::GenHit(G4int id, G4double time) :
    G4VHit(), id_(id), track_id_(0), time_(time), energy_(0), total_energy_(0), position_(0), particle_() {}

GenHit::~GenHit() {}
