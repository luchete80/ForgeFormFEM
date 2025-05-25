#include "Domain_d.h"
#include <iostream>
#include "VTKWriter.h"

#include "Mesh.h"
#include "WallTimer.h"

// #ifdef BUILD_REMESH
// #include "ReMesher.h"
// #endif

using namespace std;

namespace MetFEM{

void host_ Domain_d::Solve(){
  WallTimer timer;

  //AssignMatAddress();

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

  //CalcElemInitialVol(); //ALSO CALC VOL

  //CalcElemVol();
  //printf("calc dens\n");
  //calcElemDensity();

  //CalcNodalVol(); //To calc nodal mass
  //CalcNodalMassFromVol(); //Repla


  
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

  // /////AFTER J AND DERIVATIVES
  // if ( step_count % m_remesh_interval == 0 && step_count  >0 )
  // //if (0) //debug
  // {
    // //cout << "REMAINING " <<(step_count) % m_remesh_interval<<"INTERVAL "<<m_remesh_interval<<endl;
    // cout << "step_count "<<step_count<<endl;
    // double max=0.0;
    // int emin;
    // for (int e=0;e<m_elem_count;e++)
      // if (pl_strain[e]>max){
        // max = pl_strain[e];
        // emin = e;
      // }
    
  // //////////////////////////// IF REMESH
      // //#########################################################
      // //cout << "REMESHING "<<endl;
      // ReMesher remesh(this);
      // remesh.m_type = MMG;
      // //remesh.Generate_omegah();
      // remesh.Generate_mmg();
      // remesh.WriteDomain(); 
      // //cout << "Step "<<step_count<<endl;
      // //parallel_for ()

      // //TO MODIFY
      // double mat_cs = sqrt(mat[0]->Elastic().BulkMod()/rho[0]);
      // //SetDT(0.1*dt);
      // //dt *=0.4;

      // //double dt = 0.800e-5;
      // //cout << "New Time Step "<<dt<<endl;
      // //SetDT(dt); 
      // calcMinEdgeLength();
      // double minl = getMinLength();
      // double dt = 0.05*minl/(mat_cs);
      // cout << "min length "<<minl<<", dt "<<dt<<endl;
      // //cout << "DONE REMESH"<<endl;
      // std::string s = "out_remesh_"+std::to_string(step_count)+".vtk";
      // VTKWriter writer3(this, s.c_str());
      // writer3.writeFile();
      // remesh_ = true;

      // //#########################################################
  // //////////////////////////// IF REMESH


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
    
};
