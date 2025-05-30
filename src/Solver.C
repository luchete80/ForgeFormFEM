#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

#include "Mesh.h"
#include "WallTimer.h"

#include "Matrix.h"

#include "Solver_Eigen.h"

// #ifdef BUILD_REMESH
// #include "ReMesher.h"
// #endif

using namespace std;

namespace MetFEM{

void host_ Domain_d::Solve(){
  WallTimer timer;

  AssignMatAddress();

  cout << "done"<<endl;


  cout << "Imposing BCS"<<endl;
  
  
  // InitValues();
  
  // for (int d=0;d<m_dim;d++){
    
      // for (int n=0;n<m_node_count*m_dim;n++){
        // v[n]=a[n]=u[n]=0.0;
      // }
       // ImposeBCV(d);
   
  // }

  calcElemJAndDerivatives();
  
  //FOR PRESSURE ANP
  //CalcElemInitialVol(); //ALSO CALC VOL

  CalcElemVol();

  //IMPLICIT 
  CalcMaterialStiffElementMatrix();


  //printf("calc dens\n");
  //calcElemDensity();

  //CalcNodalVol(); //To calc nodal mass
  //CalcNodalMassFromVol(); //Repla

  double dt = 1.0e-3; 
  
  Time = 0.0;
  int step_count = 0;
  double tout = 0;
  
  bool remesh_ = false;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// MAIN SOLVER LOOP /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Main Loop----"<<endl;
  //while (Time < end_t) {
      
  if (step_count % 10000 == 0)
    printf("Step %d, Time %f\n",step_count, Time);  

  
  


    // for each element:
    // # 1. Get current nodal positions x = X + u
    // x = get_current_coords(element)     # 4 x 3
    // u = get_current_displacement(element)  # 4 x 3

    // # 2. Compute volume (Ve) and shape function gradients in current config
    // gradN, Ve = compute_gradN_and_volume(x)  # gradN: list of ∇Na (3 x 1 each)

    // # 3. Compute deformation gradient F
    // F = np.zeros((3,3))
    par_loop(e, m_elem_count){
      Matrix F(m_dim,m_dim); //DEFORMATION GRADIENT
      Matrix X(m_dim,1);
      Matrix gradN(1,m_dim);
      
      for (int n=0;n<m_nodxelem;n++){
        for (int d=0;d<m_dim;d++){
          X.Set(d,0,x[e*m_nodxelem+d]);
          gradN.Set(0,d,getDerivative(e, 0, d, n));
        }
        F += MatMul(X,gradN);
        
      }
    
    //~ // for a in range(4):
        //~ // F += np.outer(x[a], gradN[a])   # F = sum_a (x_a ⊗ ∇N_a)

    // # 4. Compute Almansi strain:
    Matrix b(3,3); //CAUCHY Green
    b = MatMul(F, F.Transpose());
    // C = F.T @ F                      # Right Cauchy-Green (not used in UL)
    // b = F @ F.T                      # Left Cauchy-Green
    // e_almansi = 0.5 * (np.eye(3) - np.linalg.inv(b))
    Matrix e(3,3);
    //e = 0.5 * (Identity(3) - b.Inv()); //Crash
    
    //CalcStressStrain(dt); /////Setting sigma
    
    // # 5. Compute stress from strain (if elastic) or plastic correction
    // sigma = constitutive_model(F, element.state_vars)   # returns Cauchy stress

    // # 6. Build B-matrix (B_mat) for material stiffness
    // B_mat = np.zeros((6, 12))   # 6 strain comp x 12 DOF (3 per node)
    // for a in range(4):
        // dN = gradN[a]  # 3x1
        // B_mat[:, 3*a:3*a+3] = [
            // [dN[0],     0,     0],
            // [0,     dN[1],     0],
            // [0,         0, dN[2]],
            // [dN[1], dN[0],     0],
            // [dN[2],     0, dN[0]],
            // [0,     dN[2], dN[1]]
        // ]

    // # 7. Compute internal force vector
    // fint = Ve * B_mat.T @ stress_to_voigt(sigma)

    // # 8. Material stiffness Kmat
    // D_tangent = compute_material_tangent(F, element.state_vars)
    // Kmat = Ve * B_mat.T @ D_tangent @ B_mat

    // # 9. Geometric stiffness Kgeo
    // Kgeo = np.zeros((12, 12))
    // for a in range(4):
        // for b in range(4):
            // Kab = gradN[a].T @ sigma @ gradN[b] * Ve
            // Kgeo[3*a:3*a+3, 3*b:3*b+3] += Kab * np.eye(3)

    // # 10. Assemble fint, Kmat, Kgeo into global system
    // assemble_global(fint, Kmat + Kgeo, element)
    
  }//ELEMENT LOOP




  cout << "Writing output "<<endl;

  // //cout << "Writing output"<<endl;
  // VTKWriter writer2(this, "out.vtk");
  // writer2.writeFile();

  // cout << "Done."<<endl;

  // ///// IF REMESH
  // /////#####################
  // //remesh.WriteDomain();
  // //calcElemJAndDerivatives();    
  
  // //////////////////////////////////////
  
  // VTKWriter writer3(this, "out_remesh.vtk");
  // writer3.writeFile();
  
  // //AFTER WRITE

  // timer.stop();
  // std::cout << "Overall elapsed time: " << timer.elapsed() << " seconds\n";  
  
  }//SOLVE
  
  
void host_ Domain_d::ElasticSolve(){
  WallTimer timer;

  AssignMatAddress();
  calcElemJAndDerivatives();
  CalcElemVol();

  //IMPLICIT 
  CalcMaterialStiffElementMatrix();

  par_loop(e, m_elem_count){
    
  }
  
  Solver_Eigen *solver = new Solver_Eigen();
  m_solver = solver;
  m_solver->setDomain(this);
  m_solver->Allocate();
  cout << "Assemblying matrix "<<endl;
  m_solver->assemblyGlobalMatrix();

  
  m_solver->applyDirichletBCs();
  cout << "Solving system"<<endl;
  m_solver->Solve();
  
  
  }//ELASTICSOLVE
    
};
