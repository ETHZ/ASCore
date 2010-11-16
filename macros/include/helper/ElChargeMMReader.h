/* 
 * File:   ElChargeMMReader.h
 * Author: hbakhshi
 *
 * Created on November 2, 2010, 2:13 PM
 */

#ifndef PMMREADER_H
#define	PMMREADER_H

#include <exception>

using namespace std;
class ElChargeMMReaderException : std::exception{
public:
    ElChargeMMReaderException(string msg) throw() : message(msg){
    };
    virtual ~ElChargeMMReaderException() throw(){};
    virtual const char* what() const throw(){
        return ("PMM READER : " + message).c_str();
    }
private:
    string message;
};

class ElChargeMMReader {
public:

    enum ChargeSource {
        GSFCharge,
        MajorityCharge,
        AllChargesMatch
    };

    enum Cuts {
        NoCUT = 0,
        D0 = 1,
        ConvRej1 = 2,
        ConvRej0 = 4,
        CTFMatch = 8
    };

    struct ValueError{
        double Value;
        double Error;

        double RelativeError() const{
            return 100.0*Error/Value;
        }
        void Print() const{
            std::cout << "Val : " << Value << " Err : " << Error << "(" << RelativeError() << "%)" << std::endl;
        }
    };
    //cuts is a bitwise or of the members of Cuts enum,
    //for example cuts can be : Cuts::D0 | Cuts::ConvRej1
    ElChargeMMReader(std::string fileName, ChargeSource cs = GSFCharge, int cuts = (D0 | ConvRej1 | CTFMatch), bool separateBarrelEndCap = true);
    virtual ~ElChargeMMReader();

    const ValueError GetPMM(double eta = 0);
private:
    TH1D hPMM;
    bool SeparateBarrelEndCap;
};

#endif	/* PMMREADER_H */

