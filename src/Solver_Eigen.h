#include "Solver.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include "Domain_d.h"
#include <unordered_map>  // For triplet storage


namespace MetFEM{
  
class Solver_Eigen:
public Solver{

public: 
  virtual int Solve();
  virtual ~Solver_Eigen(){}
  virtual void assemblyGlobalMatrix();
protected:
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::Triplet<double> T;

};





void Solver_Eigen::assemblyGlobalMatrix() {
    // Triplet format: (row, col, value)
    std::vector<std::tuple<int, int, double>> triplets;

    for (int e = 0; e < m_dom->m_elem_count; ++e) {
        
        Matrix *Ke = m_dom->m_Kmat[e];
        std::vector<int> global_dofs(m_dom->m_nodxelem * m_dim);
        for (int a = 0; a < m_dom->m_nodxelem; ++a) {
            int node = m_dom->getElemNode(e, a);
            for (int i = 0; i < m_dim; ++i) {
                global_dofs[a * m_dim + i] = node * m_dim + i;
            }
        }

        // Collect non-zero entries
        for (int i = 0; i < global_dofs.size(); ++i) {
            int I = global_dofs[i];
            for (int j = 0; j < global_dofs.size(); ++j) {
                int J = global_dofs[j];
                if (Ke->getVal(i, j) != 0.0) {  // Skip zeros
                    triplets.emplace_back(I, J, Ke->getVal(i, j) );
                }
            }
        }
    }

    // Build sparse matrix (e.g., Eigen::SparseMatrix)
    //sparse_K.setFromTriplets(triplets.begin(), triplets.end());
}


int Solver_Eigen::Solve(){

    int n = 3;
    std::vector<T> triplets;
    triplets.reserve(7);
      
      // Build the matrix A (3x3)
    triplets.push_back(T(0,0,3)); triplets.push_back(T(0,1,2)); triplets.push_back(T(0,2,-1));
    triplets.push_back(T(1,0,2)); triplets.push_back(T(1,1,-2)); triplets.push_back(T(1,2,4));
    triplets.push_back(T(2,0,-1)); triplets.push_back(T(2,1,0.5)); triplets.push_back(T(2,2,-1));

    SpMat A(n,n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd b(n);
    b << 1, -2, 0;

    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    if(solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed\n";
        return -1;
    }

    Eigen::VectorXd x = solver.solve(b);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Solving failed\n";
        return -1;
    }

    std::cout << "Solution:\n" << x << std::endl; // Should print [1, -1, -2]

    return 0;
  
  
}


//~ int main() {
    //~ typedef Eigen::SparseMatrix<double> SpMat;
    //~ typedef Eigen::Triplet<double> T;

    //~ int n = 3;
    //~ std::vector<T> triplets;
    //~ triplets.reserve(7);

    //~ // Build the matrix A (3x3)
    //~ triplets.push_back(T(0,0,3)); triplets.push_back(T(0,1,2)); triplets.push_back(T(0,2,-1));
    //~ triplets.push_back(T(1,0,2)); triplets.push_back(T(1,1,-2)); triplets.push_back(T(1,2,4));
    //~ triplets.push_back(T(2,0,-1)); triplets.push_back(T(2,1,0.5)); triplets.push_back(T(2,2,-1));

    //~ SpMat A(n,n);
    //~ A.setFromTriplets(triplets.begin(), triplets.end());

    //~ Eigen::VectorXd b(n);
    //~ b << 1, -2, 0;

    //~ Eigen::SparseLU<SpMat> solver;
    //~ solver.analyzePattern(A);
    //~ solver.factorize(A);

    //~ if(solver.info() != Eigen::Success) {
        //~ std::cerr << "Factorization failed\n";
        //~ return -1;
    //~ }

    //~ Eigen::VectorXd x = solver.solve(b);
    //~ if(solver.info() != Eigen::Success) {
        //~ std::cerr << "Solving failed\n";
        //~ return -1;
    //~ }

    //~ std::cout << "Solution:\n" << x << std::endl; // Should print [1, -1, -2]

    //~ return 0;
//~ }

};
