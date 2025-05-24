#include "Domain_d.h"
#include <iostream>
#include <vector>

#include "Matrix.h"
#include <sstream>
#include <fstream> 
#include <iostream>


#include "tensor3.C"
#include "../lib/LSDynaReader/src/lsdynaReader.h"


using namespace std;
using namespace LS_Dyna;

namespace MetFEM {
  
  

dev_t void Domain_d::calcElemJAndDerivatives () {

  par_loop (e, m_elem_count) {

  Matrix jacob(m_dim, m_dim);
  Matrix inv_j(m_dim, m_dim);
  Matrix x2(m_nodxelem, m_dim);  
  Matrix dHxy_detJ_loc(m_dim, m_nodxelem);
  

  int offset = m_gp_count * e;


   
   // //printf("Jacob\n");jacob.Print();
   double gpc[8][3];

	//printf ("Matrices created\n");

      // do i=1,nodxelem
          // !print *, "elnod " , elem%elnod(e,i)
          // x2(i,:)=nod%x(elem%elnod(e,i),:)
      // end do
  int nind = e * m_nodxelem;
  for (int i=0;i<m_nodxelem;i++){
      ////TEMPLATIZE
      if (m_dim == 2){
        double2 x_ = Ptr_vector2(x,m_elnod[nind+i]);
        x2.Set(i,0,x_.x); x2.Set(i,1,x_.y); 
        printf("x2\n");
        x2.Print();
      } else {
        vector_t x_ = Ptr_vector_t(x,m_elnod[nind+i]); 
        x2.Set(i,0,x_.x); x2.Set(i,1,x_.y);        
        x2.Set(i,2,x_.z);
       
      }

        
      
      ////printf ("elnod %d, %lf %lf %lf \n",m_elnod[nind+i],x[m_elnod[nind+i]].x,x[m_elnod[nind+i]].y,x[m_elnod[nind+i]].z);
  } 
  //printf("x2\n");x2.Print();
  //printf("m_gp_count %d\n",m_gp_count);
    //printf("Calculating jacobian\n");
    if (m_gp_count == 1 ) {      
      //invJ = adj(elem%jacob(e,gp,:,:)) !!! IN FACT IS invJ x detJ
			if (m_dim == 2) {
      // if (dim .eq. 2) then 
        // !dHdrs [-1,1,1,-1;  -1.-1,1,1] x X2
        // !! J = [
        // !! dx/dr dy/dr
        // !! dx/ds dy/dx ]
        // !!! THIS IS TO AVOID MATMUL
        // ! print *, "nodes X ", x2(:,1)
        // ! print *, "nodes Y ", x2(:,2)
        if (m_nodxelem == 4){
          for (int d=0;d<2;d++){
            // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)
            // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)
            // elem%jacob(e,gp,:,:) = 0.25*elem%jacob(e,gp,:,:)
            jacob.Set(0,d,0.25*(-x2.getVal(0,d)+x2.getVal(1,d)+x2.getVal(2,d)-x2.getVal(3,d))); 
            jacob.Set(1,d,0.25*(-x2.getVal(0,d)-x2.getVal(1,d)+x2.getVal(2,d)+x2.getVal(3,d)));
          }
        
          AdjMat(jacob, &inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
          //printf(" J ptr\n");
          //jacob.Print();
          //printf("ADJ J ptr\n");
          //inv_j.Print();          //printf("jacob\n");jacob.Print();
          //invj x dHdrs [-1,1,1,-1;  -1.-1,1,1] 
          for (int d=0;d<2;d++){        
            dHxy_detJ_loc.Set(d,0,0.25*(-inv_j.getVal(d,0)-inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,1,0.25*(inv_j.getVal(d,0)-inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,2,0.25*( inv_j.getVal(d,0)+inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,3,0.25*(-inv_j.getVal(d,0)+inv_j.getVal(d,1)));     
          }
        //dHxy_detJ_loc.Mul(0.25);
        } else if (m_nodxelem == 3){ //TRIANGLE CONSTANT ELEMENT
          //BENSON 2.4.5.2 N1 = r , N2 = s, n3 = 1 - r -s
           //dHdrs [1,0,-1;  -0,1,-1] x X2
          for (int d=0;d<2;d++){
            jacob.Set(0,d,(x2.getVal(0,d)-x2.getVal(2,d))); 
            jacob.Set(1,d,(x2.getVal(1,d)-x2.getVal(2,d)));          
          }
          AdjMat(jacob, &inv_j);
          //invj x dHdrs [-1, 1, 0,-1;  
          //              -1, 0, 1, 0]
          //                  
          for (int d=0;d<2;d++){    //col of dHdrs     
            dHxy_detJ_loc.Set(d,0,(inv_j.getVal(d,0)));   //row 1 of jacobian  
            dHxy_detJ_loc.Set(d,1,(inv_j.getVal(d,1)));     
            dHxy_detJ_loc.Set(d,2,(-inv_j.getVal(d,0)-inv_j.getVal(d,1)));      
          }          
        }//TRIANGLE
			} else { //!!!DIM 3
          if (m_nodxelem==8){
            for (int d=0;d<m_dim;d++){ //HEXA
              jacob.Set(0,d,0.125*(-x2.getVal(0,d)+x2.getVal(1,d)+x2.getVal(2,d)-x2.getVal(3,d)-x2.getVal(4,d)+x2.getVal(5,d)+x2.getVal(6,d)-x2.getVal(7,d)));  
              jacob.Set(1,d,0.125*(-x2.getVal(0,d)-x2.getVal(1,d)+x2.getVal(2,d)+x2.getVal(3,d)-x2.getVal(4,d)-x2.getVal(5,d)+x2.getVal(6,d)+x2.getVal(7,d)));  
              jacob.Set(2,d,0.125*(-x2.getVal(0,d)-x2.getVal(1,d)-x2.getVal(2,d)-x2.getVal(3,d)+x2.getVal(4,d)+x2.getVal(5,d)+x2.getVal(6,d)+x2.getVal(7,d))); 
              //jacob.Set(0,d,-x2.getVal(0,d) + x2.getVal(1,d) + x2.getVal(2,d) - x2.getVal(3,d));  

            }


          AdjMat(jacob, &inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
          //printf("ADJ J ptr\n");
          //inv_j.Print();          //printf("jacob\n");jacob.Print();
                  
          // jacob.Print();
          ////printf("INV J2 not ptr\n");
          //inv.Print();
          
          //inv.Print();
          
          for (int d=0;d<m_dim;d++){            
            dHxy_detJ_loc.Set(d,0,0.125*(-inv_j.getVal(d,0)-inv_j.getVal(d,1)-inv_j.getVal(d,2)));         
            dHxy_detJ_loc.Set(d,1,0.125*( inv_j.getVal(d,0)-inv_j.getVal(d,1)-inv_j.getVal(d,2)));  
            dHxy_detJ_loc.Set(d,2,0.125*( inv_j.getVal(d,0)+inv_j.getVal(d,1)-inv_j.getVal(d,2)));  
            dHxy_detJ_loc.Set(d,3,0.125*(-inv_j.getVal(d,0)+inv_j.getVal(d,1)-inv_j.getVal(d,2)));             
            dHxy_detJ_loc.Set(d,4,0.125*(-inv_j.getVal(d,0)-inv_j.getVal(d,1)+inv_j.getVal(d,2))); 
            dHxy_detJ_loc.Set(d,5,0.125*( inv_j.getVal(d,0)-inv_j.getVal(d,1)+inv_j.getVal(d,2)));
            dHxy_detJ_loc.Set(d,6,0.125*( inv_j.getVal(d,0)+inv_j.getVal(d,1)+inv_j.getVal(d,2)));
            dHxy_detJ_loc.Set(d,7,0.125*(-inv_j.getVal(d,0)+inv_j.getVal(d,1)+inv_j.getVal(d,2)));
          }
          //dHxy_detJ_loc.Mul(0.125); /////.DO NOT USE THIS!! --- ERRORS ---

          // // elem%dHxy_detJ(e,gp,:,:) = elem%dHxy_detJ(e,gp,:,:) * 0.125d0    
          } else if (m_nodxelem==4){ //TETRA

            //J =dr/dx
            for (int d=0;d<m_dim;d++){
              jacob.Set(0,d,x2.getVal(1,d)-x2.getVal(0,d) ); //d1-d4
              jacob.Set(1,d,x2.getVal(2,d)-x2.getVal(0,d) );            
              jacob.Set(2,d,x2.getVal(3,d)-x2.getVal(0,d) );      
            }
            //USE ADJ TO NOT DIVIDE BY DET
            //dHdr = dH/dr dr/dx
            AdjMat(jacob, &inv_j); //NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX
            //printf(" J ptr\n");
            //jacob.Print();
            //printf("ADJ J ptr\n");
            //inv_j.Print();          //printf("jacob\n");jacob.Print();
            //invj ((d,X) x dHdrs [-1,1,0,0;  
            //                     -1,0,1,0;
            //                     -1,0,0,1]
            for (int d=0;d<m_dim;d++){    
              /////ROWS OF INVJ
              dHxy_detJ_loc.Set(d,0,-inv_j.getVal(d,0)-inv_j.getVal(d,1)-inv_j.getVal(d,2));   
              dHxy_detJ_loc.Set(d,1, inv_j.getVal(d,0) );     
              dHxy_detJ_loc.Set(d,2, inv_j.getVal(d,1) );     
              dHxy_detJ_loc.Set(d,3, inv_j.getVal(d,2) );     
            }

          }//TETRA

          
      } // end if  !!!!DIM
      
      m_detJ[offset] = jacob.calcDet();
      //printf("det J %f allocated, offset %d\n",m_detJ[offset],offset);
      //if (m_detJ[offset]>1.0 || m_detJ[offset]<1.0e-3){
        //printf("--------------------- WARNGIN JACOBIAN\n");
        
        //if (e<10){
        ///printf("ELNODES %d %d %d %d \n",m_elnod[nind],m_elnod[nind+1],m_elnod[nind+2],m_elnod[nind+3]);
        //printf("det J %f allocated, offset %d\n",m_detJ[offset],offset);
        //}
        
      //}
      // elem%detJ(e,gp) = det(elem%jacob(e,gp,:,:))
    } else { //!!!!! GP > 1

    }// end if !!gp ==1

    ///// ALLOCATION
    for (int gp=0;gp<m_gp_count;gp++){
      //Domain_d::setDerivative(const int &e, const int &gp, const int &i, const int &j, const double &v)
      //setDerivative(e,gp,dHxy_detJ_loc
      for (int j=0;j<m_nodxelem;j++){
        int offset = e*(m_nodxelem * m_gp_count) + gp * m_nodxelem;
        ////printf ("Offset %d \n", offset);
        
          //m_dH_detJ_dx[offset + j                 ] = dHxy_detJ_loc.operator()(0,j);
          // m_dH_detJ_dx[offset + j] = dHxy_detJ_loc.getVal(0,j);
          // m_dH_detJ_dy[offset + j] = dHxy_detJ_loc.getVal(1,j); 
          // m_dH_detJ_dz[offset + j] = dHxy_detJ_loc.getVal(2,j);      
          setDerivative(e,gp,0,j,dHxy_detJ_loc.getVal(0,j));
          setDerivative(e,gp,1,j,dHxy_detJ_loc.getVal(1,j));
          if (m_dim ==3)
            setDerivative(e,gp,2,j,dHxy_detJ_loc.getVal(2,j));
          //printf("set der: z n %d %f\n",j, dHxy_detJ_loc.getVal(2,j));
          
      }
    }

  } // e < elem_colunt
  

}



};