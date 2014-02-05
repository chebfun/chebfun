function [L, f] = linop(N)
%LINOP Converts a chebop to a linop
% L = LINOP(N) converts a chebop N to a linop L if N is a linear operator. If N
% is not linear, an error message is returned.
%
% [L, F] = LINOP(N) returns also the affine part F of the linear chebop N such
% that L*u + F(x) = N.op(x,u).
%
% See also LINOP, CHEBOP/LINEARIZE, CHEBOP/ISLINEAR

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Tell CHEBOP/linearize() to stop if it detects nonlinearity.
linCheck = 1; 

% Linearize, thereby obtaining linearity information and a LINOP
[L, f, isLinear] = linearize(N, [], [], linCheck);

% We need the entire operator (including BCs) to be linear
isLinear = all(isLinear);

% TODO: Remove once linearity detection is ready.
warning('Linearity detection for chebops is not yet fully implemented');

% Throw an error is the chebop is nonlinear
if ~( isLinear )
    error('CHEBFUN:CHEBOP:linop:nonlinear',...
        'Chebop does not appear to be a linear operator.')
end