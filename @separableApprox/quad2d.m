function I = quad2d(f, a, b, c, d, varargin)
%QUAD2D  Complete definite integral of SEPARABLEAPPROX. 
%   I = QUAD2D( F ), returns the definite integral of a SEPARABLEAPPROX integrated
%   over its domain of definition.
% 
%   I = QUAD2D(F, a, b, c, d), returns the definite integral of a SEPARABLEAPPROX.
%   Integrated over the domangle [a b] x [c d].
% 
% See also INTEGRAL2, SUM2, INTEGRAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Is [a b c d]  subset of the domain of f ? 
dom = f.domain;  
if ( ( a < dom(1) ) || ( b > dom(2) ) || ( c < dom(3) ) || ( d > dom(4) ) )
    error('CHEBFUN:SEPARABLEAPPROX:quad2d:domain', ...
        'Can only integrate within the SEPARABLEAPPROX''s domain');
end

% Form a new SEPARABLEAPPROX and integrate. (This can be faster.)
I = integral2( restrict( f, [a, b, c, d] ) ); % Any extra arguments are ignored. 

end
