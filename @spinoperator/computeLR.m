function LR = computeLR(S, L, N, M, dt)
%COMPUTELR   Create a contour around each eigenvalue of the linear part of a 
%SPINOPERATOR.
%   LR = COMPUTELR(S, L, N, M, DT) outputs a NxM matrix to be used for 
%   the complex means. L is the linear part of the SPINOPERATOR S, discretized 
%   with N in points in each space dimension, M is the number of points to 
%   discretize the contour, and dt is the timestep.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the dimension DIM and the number of variables NVARS:
dim = S.dimension;
nVars = S.numVars;

% Roots of unity:
if ( isreal(L) == 1 )
    % Use only upper-half circle when eigenvalues are real:
    r = exp(1i*pi*((1:M) - .5)/M);
else
    r = exp(2i*pi*((1:M) - .5)/M);
end

% Move each root of unity around each entry of DT*L:
if ( dim == 2 || dim == 3 )
    L = L(:);
end
LR = dt*repmat(L, 1, M) + repmat(r, nVars*N^dim, 1);

end