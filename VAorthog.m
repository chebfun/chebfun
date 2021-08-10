function [Hes,R] = VAorthog(Z,n,varargin)  % Vand.+Arnoldi orthogonalization
%VAORTHOG  Vandermonde with Arnoldi orthogonalization.
%
% This code comes from Brubeck and Trefethen, "Lightning Stokes solver", arXiv 2020.
% For the mathematics, see Brubeck, Nakatsukasa, and Trefethen, "Vandermonde with
% Arnoldi", SIAM Review, 2021.
%
%  Input:   Z = column vector of sample points
%           n = degree of polynomial (>= 0)
%         Pol = cell array of vectors of poles (optional)
% Output: Hes = cell array of Hessenberg matrices (length 1+length(Pol))
%           R = matrix of basis vectors
M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end
% First orthogonalize the polynomial part
Q = ones(M,1); H = zeros(n+1,n);
for k = 1:n       
   q = Z.*Q(:,k);
   for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end 
   H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{1} = H; R = Q;
% Next orthogonalize the pole parts, if any
while ~isempty(Pol)
   pol = Pol{1}; Pol(1) = [];
   np = length(pol); H = zeros(np,np-1); Q = ones(M,1);
   for k = 1:np       
      q = Q(:,k)./(Z-pol(k));
      for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end 
      H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
   end
   Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
end
