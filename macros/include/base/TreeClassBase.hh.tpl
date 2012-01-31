//
// ASCII template for ntuple interface
//
// DO NOT MODIFY UNLESS YOU KNOW EXACTLY WHAT YOU ARE DOING
//
#ifndef TreeClassBase_h
#define TreeClassBase_h

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
    
class TreeClassBase {
public:
    TreeClassBase(const std::vector<std::string>& fileList);
    virtual ~TreeClassBase();
    virtual const bool GetAllByLabel(void);

    // Declaration of run leaf types
    <RUNLEAFDECLARATION>

    // Declaration of event leaf types
    <EVENTLEAFDECLARATION>
    
protected:
    virtual void  Init();
    Long64_t fCurRun;
    fwlite::ChainEvent *fEvent;    

};

#endif
