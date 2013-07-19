function out = sum(f)
%SUM   Definite integral of a SINGFUN on the interval [-1,1].
%   SUM(F) is the integral of F from -1 to 1.

% See also CUMSUM, DIFF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% The integral is divergent when at least one of f.exponents is smaller than or
% equal to -1.

if any(f.exponents <= -1)
    sl = sign(get(f.smoothPart, 'lval'));
    sr = sign(get(f.smoothPart, 'rval'));
    if all(f.exponents <= -1)
        if  sl == sr
            out = sl.*inf;
        else
            out = NaN;
        end
        
    elseif (f.exponents(1) <= -1 && sl == -1) || (f.exponents(2) <= -1 && sr == -1)
        out = -inf;
    else
        out = inf;
    end
    
    return
end

% When f have trivial exponents, the integral is nothing but the integral of its
% smooth part. Otherwise, we evaluate the integral by using Gauss-Jacobi points 
% and weights.

if all(f.exponents == 0)
    out = sum(f.smoothPart);
else
    n = length(f.smoothPart);
    [x, w] = jacpts(ceil(n/2)+1, f.exponents(2), f.exponents(1));    
    out = w*f.smoothPart.feval(x);
end

end