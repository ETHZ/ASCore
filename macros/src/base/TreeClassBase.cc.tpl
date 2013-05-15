//
// DO NOT MODIFY UNLESS YOU KNOW EXACTLY WHAT YOU ARE DOING
//
#include "base/TreeClassBase.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//____________________________________________________________________
TreeClassBase::TreeClassBase(const std::vector<std::string>& fileList)
{
    fEvent = new fwlite::ChainEvent(fileList);
    fCurRun = -1;
    Init();
}

//____________________________________________________________________
TreeClassBase::~TreeClassBase()
{
    if (!fEvent) return;
    delete fEvent;
}

//____________________________________________________________________
// Retrieve all branches ("getByLabel")
const bool TreeClassBase::GetAllByLabel(void) {

    const edm::EventBase* event = fEvent;
    bool result(true);

    // Initialize key variables
    Run         = 0;
    Event       = 0;
    LumiSection = 0;

    // Load run information if run has changed
    if ( fCurRun != event->id().run() ) {
      fCurRun = event->id().run();
      const edm::RunBase& run = fEvent->getRun();

      // Get all run handles and assign to members
      <GETRUNHANDLES>


    }

    // Get all event handles and assign to members
    <GETEVENTHANDLES>

    return result;
    
}


//____________________________________________________________________
// Called at init: define "tags"
void TreeClassBase::Init(void)
{

    <DEFINELABELS>

}
