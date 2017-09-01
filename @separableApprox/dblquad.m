function I = dblquad(f, a, b, c, d, varargin)
%DBLQUAD   Complete definite integral of SEPARABLEAPPROX. 
%   I = DBLQUAD(F, a, b, c, d), returns the definite integral of a SEPARABLEAPPROX over
%   the region [a, b, c, d].
% 
%   This function is a wrapper for quad2d.
%
% See also QUAD2D, INTEGRAL2, SUM2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% (Any extra arguments are ignored.)
I = quad2d( f, a, b, c, d ); 

end

