function psi = psiEval(l, C, LR, N, dim, nVars)
%PSIEVAL   Evaluate a psi-function.
%   PSI = PSIEVAL(L, C, LR, N, DIM, NVARS) evaluates the psi-function of index L 
%   and coefficient C with the contour LR, N grid points, in dimension DIM and 
%   with NVARS variables.
%
% See also EXPINTEG/PHIEVAL, EXPINTEG/PHIFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get a function handle to the phi-function of index L:
phi = expinteg.phiFun(l);

% Evaluate the psi-function with a contour integral:
psi = mean(C^l*feval(phi, C*LR), 2);

% Reshape it when nVars>1 or/and dim>1:
psi = reshape(psi, nVars*N, N^(dim>1), N^(dim>2));

end