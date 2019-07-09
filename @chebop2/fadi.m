function [UX, DX, VX] = fadi( A, B, M, N, p, q)
%FADI    Factored Alternating direction implicit method.
%X = FADI( A, B, M, N, p, q)  solves the Sylvester equation
%
%            A*X - X*B = M*N.'
%
% using the FADI method with ADI shifts p and q. It requires that the 
% righthand side is given in low-rank form, i.e., M and N where rhs = M*N.'.

[m, rho] = size( M ); 
n = size(N, 1); 
RankOfSolution = rho * numel(p);
UX = zeros(m, RankOfSolution);
VX = zeros(n, RankOfSolution);
DX = diag( kron(q-p, ones(1,rho)) );
Im = speye(m);
In = speye(n);
UX(:, 1:rho) = (A+p(1)*Im)\M;
VX(:, 1:rho) = (B+q(1)*In)\N;
for j = 1:numel(p)-1
    UX(:, j*rho+(1:rho)) = (A+q(j)*Im)*((A+p(j+1)*Im)\UX(:,(j-1)*rho+(1:rho)));
    VX(:, j*rho+(1:rho)) = (B+p(j)*In)*((B+q(j+1)*In)\VX(:,(j-1)*rho+(1:rho)));
end
end