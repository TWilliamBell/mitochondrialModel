A = readtable("../results/linearization.csv");

A = A(2:64,2:64);

A = table2array(A);
A(1:2,1:2) = 1e-3*A(1:2,1:2)

[U, S, V] = svd(A);

s = diag(S);
s