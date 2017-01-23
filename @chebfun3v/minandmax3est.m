function box = minandmax3est(f, N)
%MINANDMAX3EST   Estimates the range of a CHEBFUN3V.
%   BOX = MINANDMAX3EST(F) returns estimates for the minimum and maximum of each
%   component of the CHEBFUN3V F over its domain.  BOX is a vector of length
%   twice the number of components of F, containing the estimated minimum and
%   maximum of each component.
%
%   BOX = MINANDMAX3EST(F, N) returns estimates for the minimum and maximum of
%   each component of the CHEBFUN3V F over its domain, based on samples on
%   an N by N by N grid (N = 25 by default).
% 
% See also CHEBFUN3/MINANDMAX3EST.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

box = [];

if ( isempty(f) )
    return
end

if ( ( nargin < 2 ) || isempty(N) )
    % Default to N = 25:
    N = 25;
end

for jj = 1:f.nComponents
    box = [ box, minandmax3est(f.components{jj}, N) ];
end

end
