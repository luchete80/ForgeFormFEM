#include "Solver.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include "Domain_d.h"
#include <unordered_map>  // For triplet storage

using namespace std;

namespace MetFEM{
  
class Solver_Eigen:
public Solver{

public: 
  virtual int Solve();
  virtual ~Solver_Eigen(){}
  virtual void assemblyGlobalMatrix();
  virtual void Allocate();
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
            cout << "Element "<<e<<endl;
            int node = m_dom->getElemNode(e, a);
            //int node = m_dom->m_elnod[e*m_dom->m_nodxelem+a];
            for (int i = 0; i < m_dom->m_dim; ++i) {
                //global_dofs[a * m_dom->m_dim + i] = node * m_dom->m_dim + i;
            }
        }

        //~ for (int i = 0; i < ndof; ++i) {  // <<--- FIXED
            //~ int I = global_dofs[i];
            //~ for (int j = 0; j < ndof; ++j) {  // <<--- FIXED
                //~ int J = global_dofs[j];
                //~ double val = Ke->getVal(i, j);
                //~ if (val != 0.0) {
                    //~ triplets.emplace_back(I, J, val);
                //~ }
                //~ // Optional: keep or remove this debug print
            //~ }
        //~ }
    }

    //~ K.setFromTriplets(triplets.begin(), triplets.end());
}

int Solver_Eigen::Solve(){


    //R << 1, -2, 0;
    std::cout << "K norm: " << K.norm() << std::endl;
    cout << "Analyzing pattern"<<endl;
    solver.analyzePattern(K);
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
