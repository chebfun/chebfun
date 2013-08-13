function out = innerProduct(f, g)
%INNERPRODUCT   Compute the inner product of two SINGFUN objects.
%   INNERPRODUCT(F, G) returns the L2 inner product (on [-1,1]) of the two
%   SINGFUN objects F and G (conjugate linear in F).
%
% See also SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

if ( ~isa(f, 'singfun') || ~isa(g, 'singfun') )
    error('CHEBFUN:SINGFUN:innerProduct:input', ...
        'innerProduct() only operates on two SINGFUN objects.');
end

f.exponents = f.exponents + g.exponents;
f.smoothPart = f.smoothPart*g.smoothPart;

% Call sum in singfun:
out = sum(f);

% Force real output if the inputs are equal:
if ( isequal(f, g) )
    out = real(out);
end

end