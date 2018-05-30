function [r, pol, res, zer, z, f, w] = cleanup(z, f, w, Z, F, bpol) 
%CLEANUP   Attempts to remove spurious poles by removing support points given 
%   [R, POL, RES, ZER, Z, F, W] = CLEANUP(R, Z, F, W, BPOL) returns function
%   handle R, poles POL, residues RES, zeros ZER, support points Z, data values
%   F, and barycentric weights W that result after removing input support points
%   nearest to each of the components of BPOL.
%
% See also AAA.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
% Cleans up spurious poles given by bpol vec

n = numel(bpol);

for j = 1:n
    azp = abs(z-bpol(j));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    z(jj) = [];
    f(jj) = [];
end

% Remove support points z from sample set:
for jj = 1:length(z)
    F(Z == z(jj)) = [];
    Z(Z == z(jj)) = [];
end
m = length(z);
M = length(Z);

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(f);
C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(A, 0);
w = V(:,m);

% Remove support points corresponding to 0 weights
I = find(w == 0);
z(I) = [];
w(I) = [];
f(I) = [];

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, z, f, w);
[pol, res, zer] = prz(r, z, f, w);

end % End of CLEANUP().
