#include "helper/LeptJetStat.h"

using namespace std;

using std::vector;
using std::cout;
using std::endl;


// Constructor:

LeptJetStat::LeptJetStat(void)
{

  // Initialize the tables of configurations

   for (int j=0; j<indLeptmax; ++j){
     for (int k=0; k<indJetmax; ++k){
       configTable[j][k] = 0;
     }
   }

}


//------------------------------------------------------------------------------

void LeptJetStat::FillLeptJetStat(int indL[], int nJets, int nBJets) {

// Fills a configuration in the configurations table
//
// indL[4] = array with type of leptons, zero filled if <4 leptons
//   indL[i]: 1=e+, 2=e-, 3=mu+, 4=mu- for each lepton
// nJets = number of jets
// nBJets = number of b-jets
// ------------------------------------------------------
    
    if (nJets < nBJets) {
      cout << " WARNING: in FillLeptJetStat nJets < nBJets" << endl;
    }
    int index = IndexPack(indL);
    configTable[index][0]++;
    int indJet = 0;
    int indBJet = 0;
    if (nJets >= 0 && nJets <= 6){indJet = nJets+1;}
    else if (nJets > 6){indJet = 7;}
    if (indJet > 0) {configTable[index][indJet]++;}
    if (nBJets > 0 && nBJets <= 2){indBJet = nBJets+7;}
    else if (nBJets > 2){indBJet = 9;}
    if (indBJet > 0) {configTable[index][indBJet]++;}

}

//------------------------------------------------------------------------------

int LeptJetStat::IndexPack(int indL[])
{

  // pack the lepton indices into a single index number
  int index = 0;
  if (indL[0] == 5){
    index = indLeptmax - 1;
  } else if (indL[1] == 0){
    index = indL[0];
  } else if (indL[2] == 0){
    index = 4 + (indL[0]-1)*4 + indL[1];
  } else if (indL[3] == 0){
    index = 20 + (indL[0]-1)*16 + (indL[1]-1)*4 + indL[2];
//  } else if (indL[4] == 0){
  } else {
    index = 84 + (indL[0]-1)*64 + (indL[1]-1)*16 + (indL[2]-1)*4 
               + indL[3];
//   at present, use only up to 4 leptons
//  } else if (indL[5] == 0){
//    index = 340 + (indL[0]-1)*256 + (indL[1]-1)*64 + (indL[2]-1)*16
//               + (indL[3]-1)*4 + indL[4];
//  } else {
//    index = 1364 + (indL[0]-1)*1024 + (indL[1]-1)*256 + (indL[2]-1)*64
//               + (indL[3]-1)*16 + (indL[4]-1)*4 + indL[5];
  }
   
 return index;
}


void LeptJetStat::IndexUnpk(int index, int indL[])
{

  // unpack the single index number into lepton indices
  indL[0] = 0; indL[1] = 0; indL[2] = 0; indL[3] = 0;
  if (index >= indLeptmax-1){
    indL[0] = 5;
  } else if (index <= 4){
    indL[0] = index;
  } else if (index <= 20){
    index -= 5;
    indL[0] = index/4 + 1;
    indL[1] = index%4 + 1;
  } else if (index <= 84){
    index -= 21;
    indL[0] = index/16 + 1;
    indL[1] = (index%16)/4 + 1;
    indL[2] = index%4 + 1;
  } else if (index <= 340){
    index -= 85;
    indL[0] = index/64 + 1;
    indL[1] = (index%64)/16 + 1;
    indL[2] = (index%16)/4 + 1;
    indL[3] = index%4 + 1;
//   at present, use only up to 4 leptons
//  } else if (index <= 1364){
//    index -= 341;
//    indL[0] = index/256 + 1;
//    indL[1] = (index%256)/64 + 1;
//    indL[2] = (index%64)/16 + 1;
//    indL[3] = (index%16)/4 + 1;
//    indL[4] = index%4 + 1;
//  } else if (index <= 5460){
//    index -= 1365;
//    indL[0] = index/1024 + 1;
//    indL[1] = (index%1024)/256 + 1;
//    indL[2] = (index%256)/64 + 1;
//    indL[3] = (index%64)/16 + 1;
//    indL[4] = (index%16)/4 + 1;
//    indL[5] = index%4 + 1;
  } else {
    cout << " Warning *** wrong index value = " << index << endl;
  }
   
 return;
}

//------------------------------------------------------------------------------

void LeptJetStat::FillShortTables()
{

  // Fills the lepton/jet configurations in compact form from the full table
  // and prepares the correspondance list in the array indconf
  
    for (int j = 0; j < indLeptmax; ++j){
      for (int k = 0; k < indJetmax; ++k){
        shortconfig[j][k] = 0;
      }
    }
  
  // Fill the short table with the numbers of events in each configuration
    for (int j = 0; j < indLeptmax; ++j){
      int indexconf = GetConfig(j);
//    cout << " index " << j << " is conf " << indexconf << endl;
      for (int k = 0; k < indJetmax; ++k){
        shortconfig[indexconf][k] += configTable[j][k];
        if (indexconf != j){
          shortconfig[j][k] = -1;
        }
      }
    }
  
  // prepare the correspondance list
//  int nconf = ListConfigs(indconf);
  ListConfigs(indconf);
  
  return;
 }

//------------------------------------------------------------------------------

void LeptJetStat::PrintConfigsFull()
{
  char buffer[200];

  // Prints the lepton/jet configurations
  cout << " Lepton/jet configurations with ordered leptons: " << endl;
  cout << endl;
//  cout << " first buffer address = " << &buffer[0] << endl;
//  cout << " last buffer address = " << &buffer[120] << endl;
  int nch = 0;
  nch += sprintf (&buffer[nch], "%s", "     ");
  nch += sprintf (&buffer[nch], "%s", "config     ");
  nch += sprintf (&buffer[nch], "%s", " total ");
  nch += sprintf (&buffer[nch], "%s", "    0j ");
  nch += sprintf (&buffer[nch], "%s", "    1j ");
  nch += sprintf (&buffer[nch], "%s", "    2j ");
  nch += sprintf (&buffer[nch], "%s", "    3j ");
  nch += sprintf (&buffer[nch], "%s", "   4j ");
  nch += sprintf (&buffer[nch], "%s", "   5j ");
  nch += sprintf (&buffer[nch], "%s", " >=6j ");
  nch += sprintf (&buffer[nch], "%s", "   1b ");
  nch += sprintf (&buffer[nch], "%s", " >=2b ");
//  nch += sprintf (&buffer[nch], "%s", "  effic ");
  nch += sprintf (&buffer[nch], "%s", "\n");
  cout << buffer << endl;
  
  nch = 0;
  for (int i = 0; i < indLeptmax; ++i){
    nch = 0;
    int indL[4] = {0, 0, 0, 0};
    nch += sprintf (&buffer[nch], "%3i ", i);
    if (i == 0){
      nch += sprintf (&buffer[nch], "%s", "no Lepton   ");
    } else if (i == indLeptmax-1){
      nch += sprintf (&buffer[nch], "%s", ">4 Leptons  ");
    } else {
      IndexUnpk(i, indL);
      for (int j = 0; j < 4; ++j){
        if (indL[j] == 1){nch += sprintf (&buffer[nch], "%s", "e+ ");}
        else if (indL[j] == 2){nch += sprintf (&buffer[nch], "%s", "e- ");}
        else if (indL[j] == 3){nch += sprintf (&buffer[nch], "%s", "mu+");}
        else if (indL[j] == 4){nch += sprintf (&buffer[nch], "%s", "mu-");}
        else {nch += sprintf (&buffer[nch], "%s", "   ");}
      }
    }
    for (int k = 0; k < indJetmax; ++k){
      nch += sprintf (&buffer[nch], "%6i ", configTable[i][k]);
      int nprint = configTable[i][k];
      if (k < 5) {
        if (nprint > 999999) {nprint = 999999;}
        nch += sprintf (&buffer[nch], "%6i ", nprint);
      } else {
        if (nprint > 99999) {nprint = 99999;}
        nch += sprintf (&buffer[nch], "%5i ", nprint);
      }
    }
//    nch += sprintf (&buffer[nch], "%8.4f ", sumepsp[i]);
    nch += sprintf (&buffer[nch], "%s", "\n");
    cout << buffer << endl;
  }
  cout << endl;
 
 return;
}

//------------------------------------------------------------------------------

void LeptJetStat::PrintConfigs()
{
  char buffer[200];

  // Prints the lepton/jet configurations in compact form
  cout << " Lepton/jet configurations: " << endl;
  cout << endl;
  int nch = 0;
  nch += sprintf (&buffer[nch], "%s", "    ");
  nch += sprintf (&buffer[nch], "%s", "config      ");
  nch += sprintf (&buffer[nch], "%s", " total ");
  nch += sprintf (&buffer[nch], "%s", "    0j ");
  nch += sprintf (&buffer[nch], "%s", "    1j ");
  nch += sprintf (&buffer[nch], "%s", "    2j ");
  nch += sprintf (&buffer[nch], "%s", "    3j ");
  nch += sprintf (&buffer[nch], "%s", "   4j ");
  nch += sprintf (&buffer[nch], "%s", "   5j ");
  nch += sprintf (&buffer[nch], "%s", " >=6j ");
  nch += sprintf (&buffer[nch], "%s", "   1b ");
  nch += sprintf (&buffer[nch], "%s", " >=2b ");
//  nch += sprintf (&buffer[nch], "%s", "  effic ");
  nch += sprintf (&buffer[nch], "%s", "\n");
  cout << buffer << endl;
  
  int nconf = ListConfigs(indconf);
//  for (int i = 0; i < nconf; ++i){
//    cout << " indconf: nconf = " << i
//         << ", iconf = " << indconf[i] << endl;
//  }
  for (int i = 0; i < nconf; ++i){
    int index = 0;
    if (i > 0){index = indconf[i];}
    nch = 0;
    int indL[4] = {0, 0, 0, 0};
    nch += sprintf (&buffer[nch], "%3i ", i);
//    nch += sprintf (&buffer[nch], "%4i %4i %4i ", i, index, indconf[i]);
    if (i == 0){
      nch += sprintf (&buffer[nch], "%s", "no Lepton   ");
    } else if (i == nconf-1){
      nch += sprintf (&buffer[nch], "%s", ">4 Leptons  ");
    } else {
      IndexUnpk(index, indL);
      for (int j = 0; j < 4; ++j){
        if (indL[j] == 1){nch += sprintf (&buffer[nch], "%s", "e+ ");}
        else if (indL[j] == 2){nch += sprintf (&buffer[nch], "%s", "e- ");}
        else if (indL[j] == 3){nch += sprintf (&buffer[nch], "%s", "mu+");}
        else if (indL[j] == 4){nch += sprintf (&buffer[nch], "%s", "mu-");}
        else {nch += sprintf (&buffer[nch], "%s", "   ");}
      }
    }
    for (int k = 0; k < indJetmax; ++k){
      int nprint = shortconfig[index][k];
      if (k < 5) {
        if (nprint > 999999) {nprint = 999999;}
        nch += sprintf (&buffer[nch], "%6i ", nprint);
      } else {
        if (nprint > 99999) {nprint = 99999;}
        nch += sprintf (&buffer[nch], "%5i ", nprint);
      }
    }
//    nch += sprintf (&buffer[nch], "%8.4f ", sumepsp[index]);
    nch += sprintf (&buffer[nch], "%s", "\n");
    cout << buffer << endl;
  }
  cout << endl;
  
 return;
}

//------------------------------------------------------------------------------
int LeptJetStat::ListConfigs(int iconf[])
{
// Generates a list of indices (in the long table) of all physical configurations
  
  const int nLmx = 4;
  int indL[nLmx] = {0, 0, 0, 0};
  int nconf = 0;
  
  iconf[0] = 0;
  nconf++;
  
  for (int nl = 0; nl < nLmx; ++nl){
    indL[0] = 1; indL[1] = 0; indL[2] = 0; indL[3] = 0;
    for (int i = 0; i < nl; ++i){
      indL[i+1] = indL[i];
    }
//    cout << indL[0] << " " << indL[1] << " " 
//         << indL[2] << " " << indL[3] << endl;
    iconf[nconf] = IndexPack(indL);
//    cout << " ListConfigs: nconf = " << nconf
//         << ", iconf = " << iconf[nconf] << endl;
    nconf++;
    
    bool down = false;
    int il = nl;
    while (il >= 0){
      if (indL[il] < nLmx){
        indL[il]++;
        if (down){
          for (int i = il; i < nl; ++i){
            indL[i+1] = indL[i];
          }
        }
//        cout << indL[0] << " " << indL[1] << " " 
//             << indL[2] << " " << indL[3] << endl;
        iconf[nconf] = IndexPack(indL);
//        cout << " ListConfigs: nconf = " << nconf
//             << ", iconf = " << iconf[nconf] << endl;
        nconf++;
        down = false;
        il = nl;
      } else{
        down = true;
        il--;
      }
//      if (il = 0) {break;}
    }
  }
  indL[0] = 5; indL[1] = 0; indL[2] = 0; indL[3] = 0;
  iconf[nconf] = IndexPack(indL);
  nconf++;

 return nconf;
}

//------------------------------------------------------------------------------
int LeptJetStat::GetConfig(int index)
{
// Converts the index in the long table to an index in the short table
  
  const int nLeptmax = 4;
  int indL[nLeptmax] = {0, 0, 0, 0};
  int indconf = 0;
  
  if (index == 0) {return indconf;}
  else if (index >= indLeptmax-1){return indLeptmax-1;}
  IndexUnpk(index, indL);
//  cout << " before " << indL[0] << " " << indL[1] << " "
//               << indL[2] << " " << indL[3];
  
// reorder the lepton types in increasing order
  bool changed = false;
  do {
    changed = false;
    for (int i = 0; i < nLeptmax; i++){
      for (int j = i+1; j < nLeptmax; j++){
        if (indL[j] > 0 && indL[i] > indL[j]){
          int il = indL[i];
          indL[i] = indL[j];
          indL[j] = il;
          changed = true;
        }
      }
    }
  } while (changed);
//  cout << " after " << indL[0] << " " << indL[1] << " "
//               << indL[2] << " " << indL[3] << endl;
  
  indconf = IndexPack(indL);
  

 return indconf;
}

