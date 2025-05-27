Matrix A(3,3);
A.Set(0,0, 3); A.Set(0,1, 2); A.Set(0,2, -1);
A.Set(1,0, 2); A.Set(1,1, -2); A.Set(1,2, 4);
A.Set(2,0, -1); A.Set(2,1, 0.5); A.Set(2,2, -1);

Matrix b(3,1);
b.Set(0,0, 1);
b.Set(1,0, -2);
b.Set(2,0, 0);

Matrix x = SolveLinearSystem(A, b);
x.Print(); // Should print the solution

// solution is [1,-1,-2]
