#include "base/TreeReader.hh"


//____________________________________________________________________
// Go to beginning and load variables
const TreeReader& TreeReader::ToBegin() {

    fwlite::ChainEvent result = fEvent->toBegin();
    LoadAll();
    return *this;

}

//____________________________________________________________________
// Increment event number and load variables
const TreeReader& TreeReader::operator++() {

    fwlite::ChainEvent result = ++(*fEvent);
    LoadAll();
    return *this;

}

//____________________________________________________________________
const Int_t TreeReader::GetEntry(const Long64_t entry)
{
    // Read contents of entry.
    if (!fEvent) return 0;
    
    Int_t result = fEvent->to(entry);
    LoadAll(); // Retrieve all the branches
    
    return result;
}

//____________________________________________________________________
// Retrieve all branches ("getByLabel")
const bool TreeReader::LoadAll(void) {

    if ( AtEnd() ) return false;

    return GetAllByLabel();
    
}
