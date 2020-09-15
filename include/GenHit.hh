#ifndef GenHit_h
#define GenHit_h

#include <G4Allocator.hh>
#include <G4ParticleDefinition.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>
#include <globals.hh>

class GenHit : public G4VHit {
public:
    GenHit(G4int, G4double);
    virtual ~GenHit();

    G4int GetID() const {return id_;}
    G4int GetTrackID() const {return track_id_;}
    G4double GetTime() const {return time_;}
    G4double GetEnergy() const {return energy_;}
    G4double GetTotalEnergy() const {return total_energy_;}
    G4ThreeVector GetPosition() const {return position_;}
    G4ParticleDefinition* GetParticle() const {return particle_;}

    void SetTrackID(G4int trackID) { track_id_ = trackID;}
    void SetTime(G4double time) { time_ = time;}
    void AddEnergy(G4double energy) { total_energy_ += energy;}
    void SetPosition(G4ThreeVector position) { position_ = position;}
    void SetParticle(G4ParticleDefinition* particle) { particle_ = particle;}

    inline void* operator new(size_t);
    inline void operator delete(void*);

private:
    G4int id_;
    G4int track_id_;
    G4double time_;
    G4double energy_;
    G4double total_energy_;
    G4ThreeVector position_;
    G4ParticleDefinition* particle_;
};

typedef G4THitsCollection<GenHit> GenHitsCollection;

extern G4ThreadLocal G4Allocator<GenHit>* GenHitAllocator;

inline void* GenHit::operator new(size_t) {
    if(!GenHitAllocator) GenHitAllocator = new G4Allocator<GenHit>;
    return (void*)GenHitAllocator->MallocSingle();
}

inline void GenHit::operator delete(void* aHit) {
    GenHitAllocator->FreeSingle((GenHit*) aHit);
}

#endif
