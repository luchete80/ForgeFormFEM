#include "Domain_d.h"
#include "Matrix.h"

using namespace std;

namespace MetFEM{
  
void Domain_d::CalcMaterialStiffElementMatrix(){
  par_loop(e, m_elem_count){
    Matrix Bmat(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?
    for (int i=0;i<m_nodxelem;i++)
      for (int d=0;d<m_dim;d++){
        printf("i")
        Bmat.Set(m_dim,m_nodxelem*i+d,getDerivative(e, 0, i, d));
      }
    if (m_dim==3){
      for (int i=0;i<m_nodxelem;i++)
        for (int d=0;d<m_dim;d++){
          int k = d+1;if (k==m_dim) k = 0;
            //printf("i %d j %d", i, k);
            //Bmat.Set(m_dim+d,m_nodxelem*i+d,getDerivative(e, 0, i, k));
            //Bmat.Set(m_dim+d,m_nodxelem*i+k,getDerivative(e, 0, i, d));
        }
    }
  }
  
}



};