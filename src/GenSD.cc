#include "GenSD.hh"
#include "G4SDManager.hh"

GenSD::GenSD(G4String name) :
    G4VSensitiveDetector(name), hits_collection_(0), hcid_(-1) {
    G4String HCname = "genCollection";
    collectionName.insert(HCname);
}

GenSD::~GenSD() {}

void GenSD::Initialize(G4HCofThisEvent* hce) {
    hits_collection_ = new GenHitsCollection(SensitiveDetectorName, collectionName[0]);
    if(hcid_ < 0) {
        hcid_ = G4SDManager::GetSDMpointer()->GetCollectionID(hits_collection_);
    }
    hce->AddHitsCollection(hcid_, hits_collection_);
}

G4bool GenSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
    G4double edep = step->GetTotalEnergyDeposit();

    if(edep/eV < .1) return true;

    G4String type = step->GetTrack()->GetDefinition()->GetParticleType();
    if(type != "nucleus" && type != "baryon") return true;

    G4StepPoint* pre_step_point = step->GetPreStepPoint();
    G4StepPoint* post_step_point = step->GetPostStepPoint();

    G4TouchableHistory* touchable = (G4TouchableHistory*)(pre_step_point->GetTouchable());
    G4int copy_no = touchable->GetVolume()->GetCopyNo();
    G4double hit_time = pre_step_point->GetGlobalTime();
    G4int track_id = step->GetTrack()->GetTrackID();

    const G4ThreeVector kPosition = pre_step_point->GetPosition();

    G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();

    G4int ix = -1;
    for(size_t i = 0; i < hits_collection_->entries(); i++) {
        if((*hits_collection_)[i]->GetID() == copy_no) {
            ix = i;
            break;
        }
    }

    if(ix > -1) {
        GenHit* hit = (*hits_collection_)[ix];
        if(hit->GetTime() > hit_time/ns) hit->SetTime(hit_time/ns);
        hit->AddEnergy(edep/MeV);
    }
    else {
        GenHit* hit = new GenHit(copy_no, hit_time);
        hit->SetTrackID(track_id);
        hit->AddEnergy(edep/MeV);
        hit->SetPosition(kPosition/mm);
        hit->SetParticle(particle);
        hits_collection_->insert(hit);
    }

    return true;
}
