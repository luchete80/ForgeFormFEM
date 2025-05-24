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
  
#define ELNOD  4   //ORIGINALLY 8
#define FACENOD 3  //ORIGINALLY 4
#define ELFAC  4   //ORIGINALLY 6
//OLD FOT HEXA, CHANGE IT

struct Face {
    int nodes[FACENOD];
    int count; // Number of occurrences of this face
};

bool dev_t areFacesEqual(const Face& f1, const Face& f2) {
    int matchCount = 0;
    for (int i = 0; i < FACENOD; i++) {
        for (int j = 0; j < FACENOD; j++) {
            if (f1.nodes[i] == f2.nodes[j]) {
                matchCount++;
                break;
            }
        }
    }
    return matchCount == FACENOD;
}
// Add a face to the face list or increment its count if already present
void dev_t addFace(Face faceList[], int& faceCount, const Face& newFace) {
    for (int i = 0; i < faceCount; i++) {
        if (areFacesEqual(faceList[i], newFace)) {
            faceList[i].count++;
            return;
        }
    }
    // Add new face
    faceList[faceCount] = newFace;
    faceList[faceCount].count = 1;
    faceCount++;
}


// Function to add all 6 faces of a hexahedron
void dev_t addTriangleFaces(Face faceList[], int& faceCount, int element[4]) {
    // Define the 6 faces of the hexahedron
    //cout << "Element nodes "<<element[0]<<", "<<element[1]<<", "<<element[2]<<", "<<element[3]<<endl;
    Face faces[ELFAC] = {
        {{element[0], element[1], element[2]}, 0}, // Front face
        {{element[0], element[1], element[3]}, 0}, // Right face
        {{element[1], element[2], element[3]}, 0}, // Back face
        {{element[2], element[0], element[3]}, 0}, // Left face
    };

    // Add each face to the face list
    for (int i = 0; i < ELFAC; i++) {
        addFace(faceList, faceCount, faces[i]);
    }
}

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
  

}//DERIVATIVES


void Domain_d::Free(){
  free_t (x);
  free_t (v);
  free_t (a);
  free_t (u);
  free_t (u_dt);
  
  free_t (prev_a);  
	//cudaMalloc((void **)&m_f, node_count * sizeof (double) * 3);
  
  free_t (m_fi); //Internal forces
  free_t (m_fe);
  
  free_t (m_mdiag);
  free_t (m_mglob); //TODO: MAKE SPARSE. DEALLOCATED AFER DIAG CALCULATION

  free_t (T);
  free_t( m_dTedt);

 
  free_t (m_H);

  free_t (m_dH_detJ_dx);
  free_t (m_dH_detJ_dy);
  free_t (m_dH_detJ_dz);


  //DOUBLE POINTERS ATTENTION
  // free_t (m_dHrs);     //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  // free_t (x2);         //LOW ACCESS SPEED; BUT NOT DYNAMIC CREATION
  // free_t (dH_dxyz); 
  
  free_t(elnod_h); ////// USED TO COMPUTE GLOBAL M MATRIX WHICH IS COMPUTED IN CPU (TODO: CHANGE)       
  free_t(dHxy_detJ ); ////NOT USE DIRECTLY VOLUME SINCE STRAINS ARE CALC WITH THIS MATRIX

  

  free_t(m_detJ );    
  
  free_t(m_nodxelem_e);
  
  free_t(m_ematm); //Elemental mas matrices

  free_t(m_str_rate);   
  free_t(m_rot_rate );     
  free_t(m_sigma);   
  free_t(m_tau );   
  

  free_t(p); 
  free_t(rho);   
  free_t(rho_0); 
  
  free_t(vol); 
  free_t(vol_0); 

  free_t(m_voln); 
  free_t(m_vol_0n); 

  free_t(m_f_elem);   
  free_t(m_f_elem_hg);   

  free_t(mat);   

  // AXISYMM
  free_t (m_radius);

  free_t(m_jacob);
    

  free_t(ext_nodes);
  free_t(contforce);

  free_t (pl_strain);
  free_t (sigma_y);  

  free_t(m_nodel);           // NODE ELEMENTS: ELEMENTS SHARED BY EACH NODE [nodxelem* node_count] call:  m_nodel[nod,elcount] =EL i,e m_nodel[nod_offset+elcount]
  free_t(m_nodel_loc);        //
  free_t(m_nodel_offset);    //OFFSET OF THE
  free_t(m_nodel_count);    
  free_t(m_elnod_count);   /// FOR CONTACT, TO REPLACE FOR m_node_count
  
  // unsigned int    *m_contsurf_count;
  // unsigned int    *m_contsurf_elemcount;   //FOR EACH OF THE ABOVE  
  // unsigned int    *m_contsurf_elem;        //ELEMENT POS OF THE CONTACT ELEMENT 

  // ////////////////////// CONTACT 
	// // TODO, EACH RIGID PARTICLE SHOULD 
  // int   *contelem; //ELEMENT OF TRIMESH FROM "RIGID" PARTICLE, ALL FIRST PARTICLES ARE ZERO
  // TriMesh_d *trimesh;
  // int trimesh_count;
  // int *mesh_id; //particle mesh ID	
	// bool *ext_nodes;
  // int ext_nodes_count;
  // double *contforce; 
  // bool contact;
  
  //free_t(bcx_nod);  free_t(bcy_nod);  free_t(bcz_nod);
  //free_t(bcx_val);  free_t(bcy_val);  free_t(bcz_val);
  //bc_count[0] = bc_count[1] = bc_count[2] = 0;
}


void Domain_d::AssignMaterial (Material_ *material_h) {
//    cudaMalloc((void**)&materials, 1 * sizeof(Material_ )); //
    malloc_t(materials, Material_,1);
    memcpy_t(materials, material_h, 1 * sizeof(Material_));	
}

void Domain_d::setDensity(const double &r){
  double *rho_h = new double [m_elem_count];
  for (int n=0;n<m_elem_count;n++) rho_h[n] = r;
  memcpy_t(this->rho_0, rho_h, sizeof(double) * m_elem_count);    
  
  delete rho_h;
  
}

dev_t void Domain_d::SearchExtNodes() {

    printf("Adding hexas\n");
    // Array to store all faces
    //Face faceList[MAX_FACES];
    Face *faceList;
    malloc_t(faceList, Face, m_elem_count*ELFAC);
    int faceCount = 0;
    int elements[ELNOD]; //ORIGINALLY FOR HEXA IT WAS 8

    // Process each hexahedron to extract its faces
    for (int i = 0; i < m_elem_count; i++) {
        for (int ne=0;ne<m_nodxelem;ne++)
          elements[ne] = m_elnod[m_nodxelem*i+ne]; //CHANGE IF MIXED 
        //cout << "Adding faces "<<endl;
        addTriangleFaces(faceList, faceCount, elements);
    }
    //cout << "done. Face count: "<<faceCount<<endl;
    // Array to track external nodes
    for (int n=0;n<m_node_count;n++)
      ext_nodes[n] = false;
    ext_nodes_count = 0;
    //bool externalNodes[m_node_count] = {false};

    // Identify external nodes by checking faces that appear only once
    int ext_faces = 0;
    for (int i = 0; i < faceCount; i++) {
        if (faceList[i].count == 1) { // External face
            for (int j = 0; j < FACENOD; j++) {
                ext_nodes[faceList[i].nodes[j]] = true;
                
            }
            ext_faces++;
        }
    }

    // Output the external nodes
    printf("External Nodes: ");
    for (int i = 0; i < m_node_count; i++) {
        if (ext_nodes[i]) {
            printf("%d ", i);
            ext_nodes_count++;
        }
    }
    printf("\n");
    printf("Ext node count %d\n\n",ext_nodes_count);
    printf("Ext face count %d\n\n",ext_faces);
}


};