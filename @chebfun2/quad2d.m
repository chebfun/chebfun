function I = quad2d(f, a, b, c, d, varargin)
% QUAD2D Complete definite integral of chebfun2. 
%
% I = QUAD2D(F), returns the definite integral of a chebfun2. Integrated
% over its domain of definition. 
% 
% I = QUAD2D(F,a,b,c,d), returns the definite integral of a chebfun2.
% Integrated over the domangle [a b] x [c d].
% 
% See also INTEGRAL2, SUM2, INTEGRAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Is [a b c d]  subset of the domain of f ? 
dom = f.domain;  
if ( ( a < dom(1) ) || ( b > dom(2) ) || ( c < dom(3) ) || ( d > dom(4) ) )
    error('CHEBFUN2:QUAD2D','Can only integrate within the chebfun2''s domain');
end

% Be lazy and form a new chebfun2.Can be made faster. 
I = integral2( restrict(f, [a,b,c,d]) ); 

% Any extra arguments are ignored. 

end