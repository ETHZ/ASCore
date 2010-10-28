//
// Class to count selected events
//
// Example usage:
// 1) Add a few counters in class definition (with enum)
//
// #include "helper/Monitor.hh"
// ...
// private:
//   enum counters_t { count_begin, EV=count_begin, MU, EL, JET, count_end };
//   Monitor counters[count_end];
//
// 
// 2) Fill them in the analysis code
//
//    for ( ... loop over muons ... ) {
//      counters[MU].fill("All muons");
//      if ( !fTR->MuIsGMPT[index] )        return false;
//      counters[MU].fill("... is global muon prompt tight");
//      if ( !fTR->MuIsGlobalMuon[index] )  return false;
//      counters[MU].fill("... is global muon");
//      if ( !fTR->MuIsTrackerMuon[index] ) return false;
//      counters[MU].fill("... is tracker muon");
//      ...
//    }
//
//
// 3) Print out at the end
//    
//    UserAnalysis::End() {
//       ...
//       std::cout << setfill('=') << std::setw(70) << "" << std::endl;
//       std::cout << "Statistics" << std::endl;
//       std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
//       for ( counters_t iCount=count_begin; iCount<count_end; 
///            iCount = counters_t(iCount+1) ) {
//         counters[iCount].print();
//         std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
//       }
//     ...
//    }
//


#include <iomanip>
#include <string>

struct eqstr 
{
  bool operator()(const std::string& s1, const std::string s2) const
  {
    return s1.compare(s2)<0;
  }
};


class Monitor {

  typedef map<const std::string, int, eqstr> Cmap;

public:
  Monitor() : maxLength(50) {}

  void fill( const std::string& counter ) {
    Cmap::iterator it;
    if ( ( it = counters.find(counter)) == counters.end() ) {
      countNames.push_back( counter ); // Store name in ordered list
      if ( counter.length()>maxLength ) maxLength = counter.length();
      counters.insert( std::make_pair(counter,1) ); // Increment counter
    } else {
      (*it).second++;
    }
  }
  void print() {
    // This needs to be improved
    for ( std::vector<std::string>::const_iterator it = countNames.begin(); 
          it != countNames.end(); ++it ) {
      std::cout << setw(maxLength+5) << left << (*it) 
                << counters[*it] << std::endl;
    }
  }

private:
  std::vector<std::string> countNames;
  Cmap counters;
  int maxLength;

};
