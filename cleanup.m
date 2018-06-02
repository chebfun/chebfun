function [r, pol, res, zer, zj, fj, wj] = cleanup(zj, fj, wj, Z, F, bpol) 
%CLEANUP   Attempts to remove spurious poles.
%   [R, POL, RES, ZER, ZJ, FJ, WJ] = CLEANUP(ZJ, FJ, WJ, Z, F, BPOL) returns function
%   handle R, poles POL, residues RES, zeros ZER, support points ZJ, data values
%   FJ, and barycentric weights WJ that result after removing support points
%   nearest to each entry of the bad pole vector BPOL.
%
% See also AAA.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
% Cleans up spurious poles given by bpol vec

n = numel(bpol);

for j = 1:n
    azp = abs(zj-bpol(j));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    zj(jj) = [];
    fj(jj) = [];
end

% Remove support points z from sample set:
for jj = 1:length(zj)
    F(Z == zj(jj)) = [];
    Z(Z == zj(jj)) = [];
end
m = length(zj);
M = length(Z);

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(fj);
C = 1./bsxfun(@minus, Z, zj.');      % Cauchy matrix.
A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(A, 0);
wj = V(:,m);

% Remove support points corresponding to 0 weights
I = find(wj == 0);
zj(I) = [];
wj(I) = [];
fj(I) = [];

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, zj, fj, wj);
[pol, res, zer] = prz(r, zj, fj, wj);

end % End of CLEANUP().
