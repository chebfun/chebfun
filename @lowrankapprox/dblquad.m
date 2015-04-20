function I = dblquad(f, a, b, c, d, varargin)
%DBLQUAD   Complete definite integral of CHEBFUN2. 
%   I = DBLQUAD(F, a, b, c, d), returns the definite integral of a CHEBFUN2 over
%   the region [a, b, c, d].
% 
%   This function is a wrapper for quad2d.
%
% See also QUAD2D, INTEGRAL2, SUM2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% (Any extra arguments are ignored.)
I = quad2d( f, a, b, c, d ); 

end

