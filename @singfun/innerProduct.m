function out = innerProduct(f, g)
%INNERPRODUCT   Inner product of two SINGFUN objects.
%   INNERPRODUCT(F, G) returns the L2 inner product (on [-1,1]) of the two
%   SINGFUN objects F and G (conjugate linear in F).
%
% See also SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Deal with empty case:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

if ( ~isa(f, 'singfun') || ~isa(g, 'singfun') )
    error('CHEBFUN:SINGFUN:innerProduct:input', ...
        'innerProduct() only operates on two SINGFUN objects.');
end

% Add the exponents:
f.exponents = f.exponents + g.exponents;
% Multiply the smoothparts:
f.smoothPart = f.smoothPart.*conj(g.smoothPart);

% Call SUM in SINGFUN:
out = sum(f);

end
