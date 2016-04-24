function LR = computeLR(S, dt, L, M)
%COMPUTELR   Create a contour around each eigenvalue of the linear part of a 
%SPINOPERATOR.
%   LR = COMPUTELR(S, DT, L, M) outputs a matrix to be used for the complex
%   means. DT is the time-step, L is the linear part of the SPINOPERATOR S, and
%   M is the number of points to discretize the contour.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dim = getDimension(S); % spatial dimension (1, 2 or 3)
nVars = S.numVars;     % number of unknown functions
N = size(L, 1)/nVars;  % grid points

% Roots of unity:
if ( isreal(L) == 1 )
    % Use only upper-half circle when eigenvalues are real:
    r = exp(1i*pi*((1:M) - .5)/M);
else
    r = exp(2i*pi*((1:M) - .5)/M);
end

% Move each root of unity around each entry of DT*L:
LR = dt*repmat(L(:), 1, M) + repmat(r, nVars*N^dim, 1);

end