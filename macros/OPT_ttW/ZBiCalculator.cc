#include "ZBiCalculator.h"

#include "TMath.h"


float ZBiCalculator::computeZBi( float obs, float b_pred, float b_pred_err ) {

  float tau = b_pred / ( b_pred_err*b_pred_err );
  float n_off = tau*b_pred;
  float P_Bi = TMath::BetaIncomplete( 1./(1.+tau), obs, n_off+1. );
  float Z_Bi = sqrt(2)*TMath::ErfInverse( 1 - 2.*P_Bi );

  return Z_Bi;

}
