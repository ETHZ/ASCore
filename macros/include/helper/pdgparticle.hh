#ifndef pdgparticle_hh
#define pdgparticle_hh

#include <stdlib.h>
#include <iostream>

using namespace std;
class pdgparticle{

	int		fPDGID;	
	string	fName;
	string	fTexName;
	
	int		fType;
	string 	fTypeName;
	
	public:
		pdgparticle(){}
		pdgparticle(int pdgid, string Name, string Texname){
			fPDGID = pdgid;
			fName = Name;
			fTexName = Texname;
		};
		pdgparticle(int pdgid, string Name, string Texname, int type, string Typename){
			fPDGID = pdgid;
			fName = Name;
			fTexName = Texname;
			fType = type;
			fTypeName = Typename;
		};
		
		void print(){
			cout << "-----------------------------------" << endl;
			cout << "PDG ID   : " << fPDGID << endl;
			cout << "Name     : " << fName << endl;
			cout << "TeXName  : " << fTexName << endl;
			cout << "Type     : " << fType << endl;
			cout << "TypeName : " << fTypeName << endl;

		}

		void edit_pdgid(int newpdgid){
			fPDGID = newpdgid;
		}

		void edit_name(string newname){
			fName = newname;
		}

		void edit_texname(std::string newtexname){
			fTexName = newtexname;
		}

		void edit_type(int newtype){
			fType = newtype;
		}

		void edit_typename(std::string newtypename){
			fTypeName = newtypename;
		}
		
		int get_pdgid(){
			return fPDGID;
		}
		
		string get_name(){
			return fName;
		}

		string get_texname(){
			return fTexName;
		}

		int get_type(){
			return fType;
		}

		string get_typename(){
			return fTypeName;
		}

	private:
};

#endif
