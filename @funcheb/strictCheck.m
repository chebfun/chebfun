function [ishappy, epslevel, cutoff] = strictCheck(f, pref)
%STRICTCHECK  Attempt to trim trailing Chebyshev coefficients in a FUNCHEB.
%   [ISHAPPY, CUTOFF, EPSLEVEL] = STRICTCHECK(F) returns an estimated location
%   the CUTOFF at which the FUNCHEB F could be truncated to maintain an accuracy
%   of EPSLEVEL relative to F.vscale and F.hscale. ISHAPPY is logical TRUE if
%   CUTOFF < min(length(F.values),2) or F.vscale = 0, and FALSE otherwise.
%
%   [ISHAPPY, CUTOFF, EPSLEVEL] = STRICTCHECK(F, PREF) allows additional
%   preferences to be passed. In particular, one can adjust the target accuracy
%   with PREF.funcheb.EPS.
%
%   STRICTCHECK tests to see if the absolute values of the entries in the tail
%   of the coeffs, i.e., f.coeffs(1:TESTLENGTH,:), where
%       TESTLENGTH = n,             for n = 1:3
%       TESTLENGTH = 3,             for n = 4:25
%       TESTLENGTH = round((n-1)/8) for n > 25
%   all lie below the value in PREF.funcheb.eps. This value is returned in
%   EPSLEVEL and CUTOFF is the location of the first entry above EPSLEVEL in
%   absolute value.
%
% See also strictCheck.m, looseCheck.m.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with special cases -------------------------------------------------------

% Determine n (the length of the input).
n = length(f);

% Assume we're not happy. (N'aww! :( )
ishappy = false; 

% Grab some preferences:
if ( nargin == 1 )
    pref = f.pref();
    epslevel = pref.funcheb.eps;
elseif ( isnumeric(f) )
    epslevel = pref;
%     pref = f.pref();
else
    epslevel = pref.funcheb.eps;
end

% Deal with the trivial case 
if ( n < 2 ) % (Can't be simpler than a constant!)
    cutoff = n;
    return
end

% Check the vertical scale:
if ( max(f.vscale) == 0 )
    % This is the zero function, so we must be happy!
    ishappy = true;
    cutoff = 1;
    return
elseif ( any(isinf(f.vscale)) )
    % Inf located. No cutoff.
    cutoff = n;
    return
end

% NaNs are not allowed
if ( any(isnan(f.coeffs(:))) )
    error('CHEBFUN:FUN:strictCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

% Check for convergence and chop location --------------------------------------
testLength = min(n, max(5, round((n-1)/8))); 
ac = bsxfun(@rdivide, abs(f.coeffs), f.vscale);
f.coeffs(ac < epslevel) = 0;
tail = f.coeffs(1:testLength,:);
if ( ~any(tail(:)) )
    cutoff = n - find(max(f.coeffs, [], 2) > 0, 1, 'first') + 1;
    ishappy = true;
    return
end
