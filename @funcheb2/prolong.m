function g = prolong(g, nOut)
%PROLONG    Manually adjust the number of points used in a FUNCHEB2.
%
%   GOUT = PROLONG(G, NOUT) returns a FUNCHEB2 GOUT where length(GOUT) = NOUT.
%   GOUT represents the same function as G, but using more interpolation
%   points/Chebyshev coefficients then were stored in G.
%
%   If NOUT < length(G) than the representation is compressed (by aliasing the
%   coefficients), which may result in loss of accuracy.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Store the number of values the input function has:
nIn = size(g.values, 1);

% m is the number of new values needed (negative if compressing).
nOutMinusnIn = nOut - nIn;

% Trivial case
if ( nOutMinusnIn == 0 )
    % Nothing to do here!
    return
end

% Constant case:
if ( nIn == 1 )
    g.values = g.values*ones(nOut, 1);
    g.coeffs = [zeros(nOut-1,1) ; g.coeffs(1)];
    return
end

% Prolong the points; Barycentric formula when n is small, FFT when large.
if ( nOutMinusnIn < 0 && nOut < 33 && nIn < 1000)
    % Use bary to compress:
    
    g.values = funcheb2.bary(funcheb2.chebpts(nOut), g.values);
    g.coeffs = funcheb2.chebpoly(g.values);
    
else
    % Use FFTs: (in alias.m)

    g.coeffs = funcheb2.alias(g.coeffs, nOut);
    g.values = funcheb2.chebpolyval(g.coeffs); 
        
end
