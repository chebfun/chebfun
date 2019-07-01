function g = gammaEval(j, k, LR, N, dim, nVars)
%GAMMAEVAL   Evaluate a gamma-function.
%   g = GAMMAEVAL(L, LR, N, DIM, NVARS) evaluates the gamma-function (J, K) with
%   the contour LR, N grid points, in dimension DIM and with NVARS variables.
%
% See also EXPINTEG/GAMMAFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get a function handle to the phi function of index L:
g = expinteg.gammaFun(j, k);

% Evaluate it with a contour integral:
g = mean(feval(g, LR), 2);

% Reshape it when nVars>1 or/and dim>1:
g = reshape(g, nVars*N, N^(dim>1), N^(dim>2));

end