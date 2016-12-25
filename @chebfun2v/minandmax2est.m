function box = minandmax2est(f, N)
%MINANDMAX2EST   Estimates the range of a CHEBFUN2V.
%   BOX = MINANDMAX2EST(F) returns estimates for the minimum and maximum of each
%   component of the CHEBFUN2V F over its domain.  BOX is a vector of length
%   twice the number of components of F, containing the estimated minimum and
%   maximum of each component.
%
%   BOX = MINANDMAX2EST(F, N) returns estimates for the minimum and maximum of
%   each component of the CHEBFUN2V F over its domain, based on samples on
%   an N by N grid (N = 33 by default).
% 
% See also CHEBFUN2/MINANDMAX2EST.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

box = [];

if ( isempty(f) )
    return
end

if ( ( nargin < 2 ) || isempty(N) )
    % Default to N = 33:
    N = 33;
end

for jj = 1:f.nComponents
    box = [ box, minandmax2est(f.components{jj}, N) ];
end

end
