function C = coeffs3(f, m, n, p)
% COEFFS3   Trivariate Cheybshev expansion coefficients of F. 
%   C = COEFFS3(F) returns the tensor of trivariate coefficients.
%
%   [core, C, R, T] = COEFFS3(f) returns the same coefficients kept in
%   the Tucker form.
%
%   X = COEFFS3(f, M, N, P) returns coefficients as an M x N x P tensor.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[mf, nf, pf] = size(f);
F = f.coeffs;

G = zeros(m,n,pf);

for k = 1:pf
    G(:,:,k) = chebtech2.alias(trigtech.alias(F(:,:,k).',n).',m);
end

G = permute(G,[3,2,1]);
C = zeros(p,n,m);

for k = 1:m
   C(:,:,k) = trigtech.alias(G(:,:,k),p); 
end

C = permute(C,[3,2,1]);
end