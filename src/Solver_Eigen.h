#include "Solver.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include "Domain_d.h"
#include <unordered_map>  // For triplet storage

#include "defs.h"

using namespace std;

namespace MetFEM{
  
class Solver_Eigen:
public Solver{

public: 
  virtual int Solve();
  virtual ~Solver_Eigen(){}
  virtual void assemblyGlobalMatrix();
  virtual void Allocate();
  virtual void applyDirichletBCs();

    //~ void setSolverDirect() { 
        //~ solver.reset(new Eigen::SparseLU<SpMat>());
    //~ }
    
    //~ void setSolverIterative(double tol = 1e-8) {
        //~ auto iterative = new Eigen::ConjugateGradient<SpMat>();
        //~ iterative->setTolerance(tol);
        //~ solver = std::move(iterative);
    //~ }
    
protected:
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::Triplet<double> T;

  
  SpMat K;
  Eigen::VectorXd R;
  
  Eigen::VectorXd U;
  
  Eigen::SparseLU<SpMat> solver;
  
};


void Solver_Eigen::Allocate(){
  m_dof = m_dom->m_node_count * m_dom->m_dim;
  cout << "Allocate for Domain DOFs: "<< m_dof<<endl;
  K.resize(m_dof,m_dof);
  U.resize(m_dof);
  R.resize(m_dof);

  // Now K has positive dimensions and you can safely call .sum(), .norm(), etc.
  std::cout << "K.sum() = " << K.sum() << std::endl;
  
}

void Solver_Eigen::assemblyGlobalMatrix() {
    int ndof = m_dom->m_nodxelem * m_dom->m_dim; // DOFs per element
    std::vector<T> triplets;
    triplets.reserve(m_dom->m_elem_count * ndof * ndof);

    for (int e = 0; e < m_dom->m_elem_count; ++e) {
        Matrix *Ke = m_dom->m_Kmat[e];
        std::vector<int> global_dofs(ndof);  // <<--- FIXED

        for (int a = 0; a < m_dom->m_nodxelem; ++a) {
            int node = m_dom->getElemNode(e, a);
            //int node = m_dom->m_elnod[e*m_dom->m_nodxelem+a];
            for (int i = 0; i < m_dom->m_dim; ++i) {
                global_dofs[a * m_dom->m_dim + i] = node * m_dom->m_dim + i;
            }
        }

        for (int i = 0; i < ndof; ++i) {  // <<--- FIXED
            int I = global_dofs[i];
            for (int j = 0; j < ndof; ++j) {  // <<--- FIXED
                int J = global_dofs[j];
                double val = Ke->getVal(i, j);
                if (val != 0.0) {
                    triplets.emplace_back(I, J, val);
                }
                // Optional: keep or remove this debug print
            }
        }
    }

    K.setFromTriplets(triplets.begin(), triplets.end());
}

//~ void Solver_Eigen::SetBCs(){
  
  //~ for (int dim=0;dim<m_dom->m_dim; dim++){
  //~ par_loop (n,m_dom->bc_count[dim]){
    //~ double val;
    //~ int dof;
    
    //~ if (dim == 0)       {
      //~ dof = m_dom->m_dim*bcx_nod[n]+dim 
      //~ val = m_dom->bcx_val[n];
    //~ } else if (dim == 1)  {m_dom->v[m_dom->m_dim*bcy_nod[n]+dim] = bcy_val[n];
    //~ } else if (dim == 2)  {m_dom->v[m_dom->m_dim*bcz_nod[n]+dim] = bcz_val[n];
    //~ }
  //~ }//node loop
  
//~ }

///// TODO: MODIFY WITH EIGEN MAP FOR PERFORMANCE: 
///// Eigen::Map<Eigen::VectorXi> bc_nodes(m_dom->bcx_nod, m_dom->bc_count[0]);
void Solver_Eigen::applyDirichletBCs() {
    // Phase 1: Identify all BC DOFs (parallel marking)
    std::vector<bool> is_bc_dof(m_dof, false);
    
    #pragma omp parallel for
    for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        int* nodes = nullptr;
        int count = m_dom->bc_count[dim];

        if (dim == 0) nodes = m_dom->bcx_nod;
        else if (dim == 1) nodes = m_dom->bcy_nod;
        else if (dim == 2) nodes = m_dom->bcz_nod;

        for (int i = 0; i < count; ++i) {
            int dof = nodes[i] * m_dom->m_dim + dim;
            is_bc_dof[dof] = true;  // Thread-safe for distinct DOFs
        }
    }

    // Phase 2: Apply BCs (serial but efficient)
    for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        int* nodes = nullptr;
        double* values = nullptr;
        int count = m_dom->bc_count[dim];

        if (dim == 0) {
            nodes = m_dom->bcx_nod;
            values = m_dom->bcx_val;
        } else if (dim == 1) {
            nodes = m_dom->bcy_nod;
            values = m_dom->bcy_val;
        } else if (dim == 2) {
            nodes = m_dom->bcz_nod;
            values = m_dom->bcz_val;
        }

        for (int i = 0; i < count; ++i) {
            int dof = nodes[i] * m_dom->m_dim + dim;
            double value = values[i];

            // Zero out the row (skip diagonal)
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, dof); it; ++it) {
                if (it.col() != dof) it.valueRef() = 0.0;
            }

            // Set diagonal and RHS
            K.coeffRef(dof, dof) = 1.0;
            R[dof] = value;
        }
    }

    // Optional but recommended for performance
    K.makeCompressed();

}
int Solver_Eigen::Solve(){


    //R << 1, -2, 0;
    std::cout << "K norm: " << K.norm() << std::endl;
    cout << "Analyzing pattern"<<endl;
    solver.analyzePattern(K);
    solver.setPivotThreshold(1e-6);  // Allow smaller pivots

    cout << "Done. "<<endl;
    cout << "Factorizing"<<endl;
    solver.factorize(K);
    

    if(solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed\n";
        return -1;
    }

    U = solver.solve(R);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Solving failed\n";
        return -1;
    }

    if (solver.info() != Eigen::Success) {
      std::cerr << "Factorization failed, info code: " << solver.info() << std::endl;
      return -1;
    }

    std::cout << "Solution:\n" << U << std::endl; // Should print [1, -1, -2]

    return 0;
  
  
}


};
