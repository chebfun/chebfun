function [P S Q1 Z1]=qzsplit(A, C)
%QZSPLIT a faster qz factorisation
%
% This is equivalent to standard qz, except we take account of symmetry to
% reduce the computational requirements of the QZ factorisation.

% Do the QZ by splitting the problem into two subproblems. 

A = full(A); A1 = A(1:2:end, 1:2:end); 
C = full(C); C1 = C(1:2:end, 1:2:end);
[P1 S1 Q1 Z1]=qz(A1, C1);
A2 = A(2:2:end, 2:2:end); C2 = C(2:2:end, 2:2:end);
[P2 S2 Q2 Z2]=qz(A2, C2);
[P S Q1 Z1] = reform(P1, P2, S1, S2, Q1, Q2, Z1, Z2);

end


function [P S Q Z]=reform(P1,P2,S1,S2,Q1,Q2,Z1,Z2)
% Recombine subproblems to form the QZ factorization. 

% initialise all the variables. 
hf1 = size(P1,1);
n = 2*hf1-1;
P = zeros(n);
S = zeros(n);
Q =zeros(n);
Z =zeros(n);


% push the subproblem back together
P(1:hf1, 1:hf1) = P1; 
P(hf1+1:end, hf1+1:end) = P2;

S(1:hf1, 1:hf1) = S1; 
S(hf1+1:end, hf1+1:end) = S2;

Q(1:hf1, 1:2:end) = Q1; 
Q(hf1+1:end, 2:2:end) = Q2;

Z(1:2:end, 1:hf1) = Z1; 
Z(2:2:end, hf1+1:end) = Z2;

end