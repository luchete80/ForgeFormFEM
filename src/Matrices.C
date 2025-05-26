#include "Domain_d.h"
#include "Matrix.h"

#include <iostream>


using namespace std;

namespace MetFEM{
  
void Domain_d::CalcMaterialStiffElementMatrix(){
  par_loop(e, m_elem_count){
    Matrix Bmat(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?
    for (int i=0;i<m_nodxelem;i++)
      for (int d=0;d<m_dim;d++){
        Bmat.Set(m_dim, m_dim*i+d, getDerivative(e, 0, d, m_dim*i));
      }
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
  }
  
}



};
