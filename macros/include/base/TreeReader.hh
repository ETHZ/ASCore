#ifndef TreeReader_h
#define TreeReader_h

//
// Interface to TreeClassBase
//


#include "DataFormats/FWLite/interface/Handle.h"   
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "base/TreeClassBase.hh"
    
class TreeReader : public TreeClassBase {
public:
  TreeReader(const std::vector<std::string>& fileList) : TreeClassBase(fileList) {}
  virtual ~TreeReader() {}

  // Load information
  virtual const Int_t GetEntry(const Long64_t entry);
  virtual const bool LoadAll(void);

  // Getters
  inline virtual const Long64_t GetEntries() const { return fEvent->size(); };

  // Looping interface
  virtual const TreeReader& ToBegin();
  inline virtual bool AtEnd() { return fEvent->atEnd(); }
  virtual const TreeReader &operator++();
  inline virtual const bool CheckLoading() const {
    return (Run>0)&&(LumiSection>0)&&(Event>0); 
  }

    
};

#endif
