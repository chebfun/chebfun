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

exps = f.exponents + g.exponents;

% The integral is divergent when at least one of exponents of the product is 
% smaller than or equal to -1. 

if any(exps <= -1)
    sl = sign(f.vals(1));
    sr = sign(f.vals(end));
    if all(exps <= -1)
        if  sl == sr
            out = sl.*inf;
        else
            out = NaN;
        end
        
    elseif (exps(1) <= -1 && sl == -1) || (exps(2) <= -1 && sr == -1)
        out = -inf;
    else
        out = inf;
    end
    
    return
end

% When f have trivial exponents, the integral is nothing but the integral of its
% smooth part. Otherwise, we evaluate the integral by using Gauss-Jacobi points 
% and weights.

if all(exps == 0)
    out = sum(f.smoothPart*g.smoothPart);
else
    n = length(f.smoothPart) + length(g.smoothPart);
    [x, w] = jacpts(ceil(n/2)+1, exps(2), exps(1));    
    out = w*f.smoothPart.feval(x);
end

% Force real output if the inputs are equal:
if ( isequal(f, g) )
    out = real(out);
end

end
