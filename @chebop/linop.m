function [L, f] = linop(N)
%LINOP Converts a chebop to a linop
% L = LINOP(N) converts a chebop N to a linop L if N is a linear operator.
% If N is not linear, then an error message is returned.
%
% [L F] = LINOP(N) returns also the affine part F of the linear chebop N
% such that L*u + F(x) = N.op(x,u).
%
% See also LINOP, CHEBOP/LINEARISE, CHEBOP/ISLINEAR, CHEBOP/DIFF

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% We will throw an error if the chebop is nonlinear
linCheck = 1; 

% Linearize, thus obtaining linearity information and a LINOP
[L, f, isLinear] = linearize(N, [], [], linCheck);

% We need the entire operator (including BCs) to be linear
isLinear = all(isLinear);

warning('Linearity detection for chebops is not yet fully implemented');

% Throw an error is the chebop is nonlinear
if ~( isLinear )
    error('CHEBOP:linop:nonlinear',...
        'Chebop does not appear to be a linear operator.')
end