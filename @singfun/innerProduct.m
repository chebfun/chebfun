function out = innerProduct(f, g)
%INNERPRODUCT   Inner product of two ONEFUN objects.
%   INNERPRODUCT(F, G) returns the L2 inner product (on [-1,1]) of the two
%   ONEFUN objects F and G (conjugate linear in F).
%
% See also SUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Deal with empty case:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

if ( ~isa(f, 'onefun') || ~isa(g, 'onefun') )
    error('CHEBFUN:SINGFUN:innerProduct:input', ...
        'innerProduct() only operates on two ONEFUN objects.');
end

% Call SUM:
out = sum(conj(f).*g);

end
