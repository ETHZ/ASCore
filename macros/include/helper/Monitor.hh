#ifndef DILEPTONANALYSIS_NTUPLEPRODUCER_HH
#define DILEPTONANALYSIS_NTUPLEPRODUCER_HH
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
//       }
//     ...
//    }
//


#include <iomanip>
#include <map>
#include <sstream>
#include <string>


struct eqstr 
{
  bool operator()(const std::string& s1, const std::string s2) const
  {
    return s1.compare(s2)<0;
  }
};


class Monitor {

  //typedef std::map<const std::string, float, eqstr> Cmap;
  typedef std::map<const std::string, double, eqstr> Cmap;

public:
  Monitor() : maxLength(50),name("") {}

  void setName( const std::string& s ) { name = s; }

  /*
  void fill( const std::string& counter, const float& weight=1. ) {
    Cmap::iterator it;
    if ( ( it = counters.find(counter)) == counters.end() ) {
      countNames.push_back( counter ); // Store name in ordered list
      if ( counter.length()>maxLength ) maxLength = counter.length();
      counters.insert( make_pair(counter,weight) ); // Increment counter
    } else {
      (*it).second += weight;
    }
  }
  */
  void fill( const std::string& counter, const double& weight=1. ) {
    Cmap::iterator it;
    if ( ( it = counters.find(counter)) == counters.end() ) {
      countNames.push_back( counter ); // Store name in ordered list
      if ( counter.length()>maxLength ) maxLength = counter.length();
      counters.insert( make_pair(counter,weight) ); // Increment counter
    } else {
      (*it).second += weight;
    }
  }

  const float counts( const std::string& counter ) { 
     return counters[counter];
  }

  friend ostream& operator<<( ostream& os, Monitor m ) {
    using namespace std;
    if ( !(m.countNames.size()>0) ) return os;
    // This needs to be improved
    float maxc = m.counters[m.countNames[0]];
    float prev = maxc;
    ostringstream maxclen; maxclen << maxc;
    size_t maxwidth = maxclen.str().length();
    os << left << "COUNTER> " << m.name << endl;
    for ( vector<string>::const_iterator it = m.countNames.begin();
          it != m.countNames.end(); ++it ) {
      float count = m.counters[*it];
      os << setw(m.maxLength+5) << left << (*it)
                << setw(maxwidth) << right << count << " "
                << setw(3) << right
                << static_cast<int>(maxc>0.?count/maxc*100.:0.) << "% "
                << setw(3) << right
                << static_cast<int>(prev>0.?count/prev*100.:0.) << "% "
                << endl;
      prev = count;
    }
    os << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;

    return os;
  }

  void print() { std::cout << *this; }

public:
  std::vector<std::string> countNames;
  Cmap counters;
  int maxLength;
  std::string name;

};
#endif
