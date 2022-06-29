function [S, L] = sphericalLAP(K,LapK, T, Q1, Q2, R, lambda)
%Compute smoother matrix (S) and Laplacian matrix (L)

I = eye(size(K,1));
Klamb = K + lambda*I;

C = Q2/(Q2'*Klamb*Q2)*Q2';
D = R\Q1'*(I-Klamb*C);

S = K*C+T*D;
L = LapK*C;

