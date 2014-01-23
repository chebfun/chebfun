function I = dblquad(f,a,b,c,d,varargin)
%DBLQUAD Complete definite integral of chebfun2. 
%
% I = DBLQUAD(F,a,b,c,d), returns the definite integral of a chebfun2 over
% the region [a, b, c, d]. 
% 
% This function is a wrapper for quad2d.
%
% See also QUAD2D, INTEGRAL2, SUM2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

I = quad2d( f, a, b, c, d ); % any extra arguments are ignored. 

end

