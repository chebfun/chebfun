function I = integral3(f, dom)
%INTEGRAL3  Triple integral of a CHEBFUN3 over its domain.
%   I = INTEGRAL3(F) returns the triple integral of the CHEBFUN3 object F.
%   This is the same as SUM3.
%
%   I = INTEGRAL3(F, [a b c d e g]) integrates F over the cuboid [a b] x [c
%   d] x [e, g] provided this cuboid is in the domain of F.
%
% See also CHEBFUN3/INTEGRAL, CHEBFUN3/INTEGRAL2, CHEBFUN3/SUM, 
% CHEBFUN3/SUM2 and CHEBFUN3/SUM3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    I = [];
    return
end

if ( nargin == 1 )    
    % Triple definite integral:
    I = sum3(f);
    
elseif ( nargin == 2 )
    if ( isa(dom, 'double') && ( numel(dom) == 6) )
        % Integral over restricted cuboid:
        g = restrict(f, dom);
        I = sum3(g);
    else
        error('CHEBFUN:CHEBFUN3:integral3:baddomain', ...
            'Domain should have six corners.');
    end
else
    error('CHEBFUN:CHEBFUN3:integral3:nargin', 'Too many input arguments.');
end

end