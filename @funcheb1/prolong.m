function g = prolong(g,nOut)
%PROLONG Manually adjust the number of points used in a FUNCHEB1.
%
%   GOUT = PROLONG(G, NOUT) return a FUNCHEB1 GOUT where length(GOUT) = NOUT 
%   (number of points). GOUT represents the same function as G, but using 
%   more interpolation points/Chebyshev coefficients then were stored in G.
%
%   If NOUT < length(G) then the representation is compressed (by aliasing
%   the coefficients), which may result in loss of accuracy.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

values = g.values;
coeffs = g.coeffs;

% Store the number of values the input function has
nIn = size(values, 1);

% m is the number of new values needed (negative if compressing)
nOutMinusnIn = nOut - nIn;

% Trivial case
if ( nOutMinusnIn == 0 )
    % Nothing to do here.
    return
end

% Constant case
if ( nIn == 1 )
    g.values = values*ones(nOut,1);
    g.coeffs = [zeros(nOut-1,1) ; coeffs(1)];
    return
end

% Prolong the points; Barycentric formula when n is small, FFT when large.
if ( nOutMinusnIn < 0 && nOut < 33 && nIn < 1000) % Use bary to compress
    
    g.values = funcheb1.bary(funcheb1.chebpts(nOut), values);
    g.coeffs = funcheb1.chebpoly(g.values);
    
else                                % Use FFTs
    
    if ( nOut == 1 )                % Reduce to just one coeff:
        
        e = ones(1, ceil(nIn/2)); 
        e(end-1:-2:1) = -1;
        g.values = e*coeffs(end:-2:1,:);
        g.coeffs = g.values;      
        
    else                            % Reduce to a lower degree:
        
        coeffs = funcheb1.alias(g.coeffs, nOut);
        g.values = funcheb1.chebpolyval(coeffs); 
        g.coeffs = coeffs;
        
    end
    
end
