#include "Domain_d.h"
#include "Matrix.h"

#include <iostream>


using namespace std;

namespace MetFEM{
  
void Domain_d::CalcMaterialStiffElementMatrix(){
  par_loop(e, m_elem_count){
    /// TODO: CHANGE FOR DIM = 2
    Matrix Bmat(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?
    for (int i=0;i<m_nodxelem;i++)
      for (int d=0;d<m_dim;d++)
        Bmat.Set(d, m_dim*i+d, getDerivative(e, 0, d, m_dim*i));
      
    if (m_dim==3){
      for (int i=0;i<m_nodxelem;i++)
        for (int d=0;d<m_dim;d++){
          int k = d+1;if (k==m_dim) k = 0;
            // d/dy d/dx 0   
            printf("i %d j %d der %d\n",m_dim+d,m_dim*i+d, k);
            printf("i %d j %d der %d\n",m_dim+d,m_dim*i+k, d);
            Bmat.Set(m_dim+d, m_dim*i+d, getDerivative(e, 0, k, i));
            Bmat.Set(m_dim+d, m_dim*i+k,getDerivative(e, 0, d, i));
        }
    }
    printf ("BMAT\n");
    Bmat.Print();
    Matrix BT(m_nodxelem* m_dim,2*m_dim);
    BT=Bmat.getTranspose();
    Matrix D(6,6);
    
    double G  = mat[e]->Elastic().G();
    double E  = mat[e]->Elastic().E();
    double nu = mat[e]->Elastic().Poisson();
    double f  = E/((1.0+nu)*(1.0-2.0*nu)); 
    D.Set(0,1, f*nu);                 D.Set(0,2, f*nu);
    D.Set(1,0, f*nu);                 D.Set(1,2, f*nu);
    D.Set(2,0, f*nu);D.Set(2,1, f*nu);
    
                
    MatMul(MatMul(BT,D),Bmat, m_Kmat[e]);
    printf("K ELEM\n");
    m_Kmat[e]->Print();
  }//element
  
}



};
