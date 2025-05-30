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
  m_dim = m_dom->m_node_count * m_dom->m_dim;
  cout << "Allocate for Domain DOFs: "<< m_dim<<endl;
  K.resize(m_dim,m_dim);
  U.resize(m_dim);
  R.resize(m_dim);

  // Now K has positive dimensions and you can safely call .sum(), .norm(), etc.
  std::cout << "K.sum() = " << K.sum() << std::endl;
  
}

void Solver_Eigen::assemblyGlobalMatrix() {
    int ndof = m_dom->m_nodxelem * m_dim;
    std::vector<T> triplets;
    triplets.reserve(m_dom->m_elem_count * ndof * ndof);  // Preallocate

    for (int e = 0; e < m_dom->m_elem_count; ++e) {
        Matrix *Ke = m_dom->m_Kmat[e];
        std::vector<int> global_dofs(ndof);

        for (int a = 0; a < m_dom->m_nodxelem; ++a) {
            int node = m_dom->getElemNode(e, a);
            for (int i = 0; i < m_dim; ++i) {
                global_dofs[a * m_dim + i] = node * m_dim + i;
            }
        }

        for (int i = 0; i < ndof; ++i) {
            int I = global_dofs[i];
            for (int j = 0; j < ndof; ++j) {
                int J = global_dofs[j];
                double val = Ke->getVal(i, j);
                if (val != 0.0) {
                    triplets.emplace_back(I, J, val);
                }
                else {
                  cout << "Error, val is ZERO"<<endl;
                  }
            }
        }
    }

    K.setFromTriplets(triplets.begin(), triplets.end());
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
