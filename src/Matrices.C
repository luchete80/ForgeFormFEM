#include "Domain_d.h"
#include "Matrix.h"

#include <iostream>


using namespace std;

namespace MetFEM{
  
  
//CALCULATE DEFORMATION GRADIENT F
/// F = I + Grad * dU


/// ALMANSI
//// e = 1/2(I-FT F-1)
 //////[dN1/dx    0      0 ]
     //[   0   dN1/dy    0]
     //[   0      0    dN1/dx ]    
     //[   dNdy  dN/dx    0 ]    

Matrix Domain_d::getBArrange(int &e){
  
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
}

     
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
    for (int d=0;d<3;d++) D.Set(d,d,1.0-nu);
    for (int d=3;d<6;d++) D.Set(d,d,(1.0-2.0*nu)/2.0);    
                
    MatMul(MatMul(BT,D),Bmat, m_Kmat[e]);
    //printf("K ELEM\n");
    //m_Kmat[e]->Print();
  }//element
  
}

/// KGEO - NO B APPROACH
//dev_t void Domain_d::CalcGeomStiffElementMatrix(){
   //# Insert Kab into 12x12 Kgeo_local at block [a][b]
 // Matrix Bgeo(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?  
  // # - Kgeo_local: 12x12 zero matrix (4 nodes, 3 DOFs each)

  // for a in range(4):  # loop over nodes
      // for b in range(4):
          // grad_Na = grad_N[a]  # (3x1)
          // grad_Nb = grad_N[b]  # (3x1)
          
          // # Compute 3x3 block: K_ab = Ve * (grad_Na^T * sigma * grad_Nb) * I3
          // # ACTUALLY: K_ab = Ve * outer(grad_Na, grad_Nb) : sigma (contracted)
          // # But your scalar coeff * I3 is a common APPROXIMATION (see notes below)
          // coeff = dot(grad_Na, sigma @ grad_Nb)  # scalar
          // Kab = coeff * Ve * I3  # (3x3)

          // # Insert Kab into Kgeo_local at position [3a:3a+3, 3b:3b+3]
          // Kgeo_local[3*a:3*a+3, 3*b:3*b+3] += Kab

  
  
//}

//// K GEO - B MATRIX APPROACH
dev_t void Domain_d::CalcGeomStiffElementMatrix(){
   //# Insert Kab into 12x12 Kgeo_local at block [a][b]
  Matrix Bgeo(2*m_dim, m_nodxelem* m_dim); // WITH m_dim==2?  

  Matrix Sigma(m_dim*2, m_dim *2);
  
// # Inputs:
// # - sigma: Cauchy stress tensor (3x3) [[σxx, σxy, σxz], [σxy, σyy, σyz], [σxz, σyz, σzz]]
// # - grad_N: List of shape function gradients [∇N1, ∇N2, ∇N3, ∇N4], each ∇Na is (3x1)
// # - Ve: Element volume

// # Step 1: Construct stress matrix Σ (6x6)
// Sigma = np.zeros((6, 6))
// # Fill diagonal blocks (3x3) of Σ
// Sigma[0:3, 0:3] = sigma[0, 0] * np.eye(3)  # σxx block
// Sigma[0:3, 3:6] = sigma[0, 1] * np.eye(3)  # σxy block
// Sigma[0:3, 6:9] = sigma[0, 2] * np.eye(3)  # σxz block
// # ... (repeat for σyy, σyz, σzz, ensuring symmetry)

// # Step 2: Build B_geo (6x12 matrix)
// B_geo = np.zeros((6, 12))
    
    for (int a=0;a<m_nodxelem;a++){
      double dNa[3];
      
    }
    // for a in range(4):
    // dN = grad_N[a]  # ∇Na (3x1)
    // B_geo[:, 3*a:3*a+3] = [
        // [dN[0], 0, 0],
        // [0, dN[1], 0],
        // [0, 0, dN[2]],
        // [dN[1], dN[0], 0],
        // [dN[2], 0, dN[0]],
        // [0, dN[2], dN[1]]
    // ]

// # Step 3: Compute K_geo (12x12)
// K_geo = Ve * B_geo.T @ Sigma @ B_geo


////// K GEOMETRIC
///Kij_GEO = bTi sigma bj I3x3 dV
/// FOR EACH NODE PAIR
///dNdNiT sig dNdJ I3x3

dev_t void Domain_d::assemblyForces(){

    //if ()
  //par_loop(n, m_node_count){
  for (int n=0;n<m_node_count;n++){
    for (int d=0;d<m_dim;d++)
      m_fi[n*m_dim + d] = 0.0;
      
      //printf("--------\n");    
      for (int e=0; e<m_nodel_count[n];e++) {
        int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
        int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
        int offset  = eglob * m_nodxelem * m_dim;

        for (int d=0;d<m_dim;d++){
          //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
          //if (n==9)
          //  printf("%6e %6e %6e\n",m_fi[n*m_dim],m_f_elem[n*m_dim+1],m_f_elem[n*m_dim+2]);
          m_fi[n*m_dim + d] += m_f_elem[offset + ne*m_dim + d];
        }
          if(m_thermal){
            T[n] += dt * m_dTedt[eglob*m_nodxelem+ne];
	  }
      }
      if (m_gp_count == 1 ) {  
        for (int e=0; e<m_nodel_count[n];e++) {
          int eglob   = m_nodel     [m_nodel_offset[n]+e]; //Element
          int ne      = m_nodel_loc [m_nodel_offset[n]+e]; //LOCAL ELEMENT NODE INDEX m_nodel_local
          int offset  = eglob * m_nodxelem * m_dim;
          ////printf("glob %d, loc %d \n",n,ne);
          for (int d=0;d<m_dim;d++){
            //atomicAdd(&m_f[m_elnod[n]*m_dim + d], m_f_elem[e*m_nodxelem*m_dim + n*m_dim + d]);
            m_fi[n*m_dim + d] -= m_f_elem_hg [offset + ne*m_dim + d];
          }
        }      
      }
      // printf ("force %f %f %f\n",m_fi[m_dim*n],m_fi[m_dim*n+1],m_fi[m_dim*n+2]);
    } // element


}//assemblyForcesNonLock

};
