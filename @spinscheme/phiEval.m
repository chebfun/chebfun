function phi = phiEval(l, LR, N, dim, nVars)
%PHIEVAL   Evaluate a phi function.
%   PHI = PHIEVAL(L, LR, N, DIM, NVARS) evaluates the phi function of index L 
%   with the contour LR, N grid points, in dimension DIM and with NVARS
%   variables.
%
% See also SPINSCHEME/PHIFUN, SPINSCHEME/PSIFUN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get a function handle to the phi function of index L:
phi = spinscheme.phiFun(l);

% Evaluate it with a contour integral:
phi = mean(feval(phi, LR), 2);

% Reshape it when nVars>1 or/and dim>1:
phi = reshape(phi, nVars*N, N^(dim>1), N^(dim>2));

end