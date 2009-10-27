#ifndef LeptJetStat_h
#define LeptJetStat_h


/*  \class LeptJetStat
*
*  Class where lepton and jet statistics is accumulated as a table
*
*  Authors: Luc Pape & Filip Moortgat      Date: August 2008
*                                          Updated: 
*
* Basically, the class LeptJetStat maintains a table 
*  with every lepton configuration (defined by the lepton species and charges)
*  including statistics on the number of jets and b-jets of each configuration
* The configuration is given event by event 
*  as an array with types (e+, mu+, e-, mu-)
*  with 1 entry for each lepton (max up to 4)
* This configuration is converted into a running index
*  defining the entry in the table
* During the event by event filling, a "long" table is constructed
*  without lepton reordering (e.g. e-mu != mu-e)
* If requested, a "short" (more physical) table is produced from the "long" one
*  where leptons are reordered
*  e.g. e-mu is summed with mu-e (and the entry for mu-e is zeroed)
*/
 
#include <vector>
#include <iostream>
#include <cmath>



class LeptJetStat {

public:

// Constructor:
LeptJetStat(void);

// Destructor:
virtual ~LeptJetStat(){
  
};


// Methods:

// to fill an event configuration in the table
// keeping the order in which the leptons are given (e.g. e-mu != mu-e)
// indL[4] = array with type of leptons, zero filled if <4 leptons in event
//   indL[i]: 1=e+, 2=mu+, 3=e-, 4=mu- for each lepton
// nJets = number of jets
// nBJets = number of b-jets
virtual void FillLeptJetStat(int indL[], int nJets, int NBJets);

// converts the full table to a shorter one
// e.g. e-mu is summed with mu-e
void FillShortTables(void);

// prints the full table
virtual void PrintConfigsFull(void);

// prints the shorter table
virtual void PrintConfigs(void);

// Some auxilliary functions:
// converts the array of leptons into an index in the (long) table
virtual int IndexPack(int indL[] );
// convert the index in the long table to the short one
virtual int GetConfig(int index);
// converts the order number of a configuration to the index in the short table
// iorder runs from 0 to 69
virtual int GetConfigfrOrder(int iorder) {
            return iorder<70? indconf[iorder] : -1;}

// gets the number of events for a given (short) lepton config for all jets
// (returns one entry per jet multiplicity)
// nperJets[0] = total, [1] = 1 jet, [2] = 2 jets, ..., [6] = >= 6 jets
//                      [7] = 1 b, [8] = >=2 b
virtual void NEntriesPerJetMult(int index, int nperJets[]) {
         int sum = 0;
         for (int i=1; i<indJetmax; ++i){
	   nperJets[i] = shortconfig[index][i];
           if (i < indJetmax-2) {
	     sum += shortconfig[index][i];
	   }
         }
	 nperJets[0] = sum;
         return;}

// gets the number of events for a given (short) lepton config
// (summed over jet multiplicities >= nJetsmin)
virtual int NumberOfEntries(int index, int nJetsmin) {
         if (nJetsmin < 1) {nJetsmin = 1;}
         int sum = 0;
         for (int i=nJetsmin; i<indJetmax-2; ++i){
           sum += shortconfig[index][i];
         }
         return sum;}


private:
 
  // names of parameter sets  
  
  // Define the parameters


//  User functions

   static const int indLeptmax = 342;
   static const int indJetmax = 10;
   static const int indConfmax = 71;
// Note1: indLeptmax can be computed from the number of permutations + 2
//  (+1 for no leptons and +1 for >=5 leptons)
// to each permutation we add 1 of the leptons to get the next ones
// the numbers of permutations with k leptons are thus
//  Nperm(k) = 4*Nperm(k-1) with Nperm(1) = 4
//  Nperm(1)=4, Nperm(2)=16, Nperm(3)=64, Nperm(4)=256, total=340
// Note2: indConfmax can be computed from the number of configs + 2
//  (+1 for no leptons and +1 for >=5 leptons)
// number of configs with k leptons is given by the number of combinations with repetition
//  Nconf(k) = (n+k-1)! / k! (n-1)!   where n=4 kinds of leptons
//  Nconf(1)=4, Nconf(2)=10, Nconf(3)=20, Nconf(4)=35, total=69

virtual void IndexUnpk(int, int indL[] );
virtual int ListConfigs(int iconf[]);

   int configTable[indLeptmax][indJetmax];
   int shortconfig[indLeptmax][indJetmax];
   int indconf[indConfmax];

};

#endif
