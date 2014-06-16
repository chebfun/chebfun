function out = islinear(N, u)
%ISLINEAR   Determine linearity of a CHEBOP.
%   OUT = ISLINEAR(N) determines the linearity of the CHEBOP N. In particular:
%       OUT(1) = 1 if N.OP is linear, 0 otherwise.
%       OUT(2) = 1 if N.LBC is linear, 0 otherwise.
%       OUT(3) = 1 if N.RBC is linear, 0 otherwise.
%       OUT(4) = 1 if N.BC is linear, 0 otherwise.
%
%   OUT = ISLINEAR(N, U) performs the linearization of N around the function U,
%   rather than the zero function.
%
% See also LINEARIZE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    % Allow LINEARIZE() to chose the linearization variable:
    u = [];
end

% Call LINEARIZE():
[ignored1, ignored2, out] = linearize(N, u, [], 0);

end
