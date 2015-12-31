function phit = phitEval(l, C, LR, N, dim, nVars)
%PHIEVAL   Evaluate a phit function.
%   PHI = PHIEVALT(L, LR, N, DIM, NVARS) evaluates the phit function of index L 
%   and coefficient C with the contour LR, N grid points, in dimension DIM and 
%   with NVARS variables.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get a function handle to the phi function of index L:
phi = spinscheme.phiFun(l);

% Evaluate it with a contour integral:
phit = mean(C^l*feval(phi, C*LR), 2);

% Reshape it when nVars>1 or/and dim>1:
phit = reshape(phit, nVars*N, N^(dim>1), N^(dim>2));

end