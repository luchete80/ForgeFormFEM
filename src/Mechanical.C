//CALCULATE INTERNAL FORCES
// F INT BT sigma

#include "Domain_d.h"
#include <iostream>
#include <vector>

// #include "tensor.cu"
#include "Matrix.h"

#if CUDA_BUILD
#include "tensor.cuh"
#else
//#include "Matrix.h"
#endif

#include "Tensor3.C"

using namespace std;

namespace MetFEM {

// 
// !!!!!! IT ASSUMES PRESSURE AND STRAIN RATES ARE ALREADY CALCULATED
// !!!!!! (AT t+1/2 to avoid stress at rigid rotations, see Benson 1992)
dev_t void Domain_d::CalcStressStrain(double dt){


  par_loop(e,m_elem_count){
      //printf("calculating sigma \n");
        // Jaumann rate terms
      tensor3 RotationRateT,SRT,RS;
      tensor3 RotRate;
      tensor3 StrRate;
      tensor3 ShearStress;
      tensor3 Sigma;

    //printf("calculating sigma %d\n", e);
    for (int gp=0;gp<m_gp_count;gp++){
      int offset_s = e * m_gp_count + gp;   //SCALAR
      int offset_t = offset_s * 6 ; //SYM TENSOR
      ShearStress = FromFlatSym(m_tau,          offset_t );
      StrRate     = FromFlatSym(m_str_rate,     offset_t );
      RotRate     = FromFlatAntiSym(m_rot_rate, offset_t );


      SRT = ShearStress * Trans(RotRate);
      RS  = RotRate * ShearStress;

      tensor3 test = StrRate-1.0/3.0*(StrRate.xx+StrRate.yy+StrRate.zz)*Identity();
      ShearStress	= ShearStress  + dt*(2.0* mat[e]->Elastic().G()*(StrRate - 1.0/3.0*Trace(StrRate) * Identity() ) + SRT+RS);
      
      double J2 = 0.5*(ShearStress.xx*ShearStress.xx +  2.0*ShearStress.xy*ShearStress.xy + 
                                      2.0*ShearStress.xz*ShearStress.xz + 
                     ShearStress.yy*ShearStress.yy+  
                                      2.0*ShearStress.yz*ShearStress.yz +               
                     ShearStress.zz*ShearStress.zz                 
                                     
                    );
      double sig_trial = sqrt(3.0*J2);

      if (sigma_y[e]<sig_trial){

        ShearStress = ShearStress * (sigma_y[e] / sig_trial);
        pl_strain[e] += (sig_trial - sigma_y[e]) / (3.0 *  mat[e]->Elastic().G());


      }

      Sigma = -p[offset_s] * Identity() + ShearStress;
 
      double Ep = 0;
			double dep=( sig_trial - sigma_y[e])/ (3.*mat[e]->Elastic().G() + Ep);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0

      ///// OUTPUT TO Flatten arrays
      ToFlatSymPtr(Sigma, m_sigma,offset_t);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
      //ToFlatSymPtr(Strain, 	strain,6*i);		
      ToFlatSymPtr(ShearStress, m_tau, offset_t);
      
    }//gp
  }//el < elcount


}

//To calculate before elemet Jacobian calc
dev_t void Domain_d::CalcElemVol(){
  par_loop(e,m_elem_count){
    double w;
    //TODO: CHANGE WEIGHT TO ARRAY
    if (m_gp_count == 1) {
      if (m_dim == 2) w = 4;//w = pow(2.0, m_dim);
      if (m_dim == 3)     
        if      (m_nodxelem == 4)  w = 1.0/6.0;
        else if (m_nodxelem == 8)  w = 8.0;
    } else                  w = 1.0;
    
    int offset = m_gp_count * e;
    vol[e] = 0.0;
    for (int gp=0;gp<m_gp_count;gp++){
      vol[e] += m_detJ[offset] * w;
    }  
    //if (e<10)
    //printf("Element %d Vol %f, det %f\n",e,vol[e],m_detJ[offset]);  
  }//el
}



};
