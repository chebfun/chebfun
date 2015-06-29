function [ishappy, epslevel, cutoff] = strictCheck(f, values, vscl, pref)
%STRICTCHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = STRICTCHECK(F, VALUES, VSCL) returns an
%   estimated location CUTOFF at which the CHEBTECH F could be truncated to
%   maintain an accuracy of EPSLEVEL relative to VSCL and F.HSCALE. ISHAPPY is
%   TRUE if CUTOFF < MIN(LENGTH(F.COEFFS), 2) or F.VSCALE=0, and FALSE
%   otherwise. If ISHAPPY is false, EPSLEVEL returns an estimate of the accuracy
%   achieved.
%
%   [ISHAPPY, EPSLEVEL, CUTOFF] = STRICTCHECK(F, VALUES, VSCL, PREF) allows
%   additional preferences to be passed. In particular, one can adjust the
%   target accuracy with PREF.EPS. The VALUES field is ignored, but included for
%   consistency with other happiness checks.
%
%   STRICTCHECK tests to see if the absolute values of the entries in the tail
%   of coeffs, i.e., f.coeffs(1:TESTLENGTH,:), where
%       TESTLENGTH = n,             for n = 1:4
%       TESTLENGTH = 5,             for n = 5:44
%       TESTLENGTH = round((n-1)/8) for n > 44
%   all lie below the value in PREF.EPS. This value is returned in
%   EPSLEVEL and CUTOFF is the location of the first entry above EPSLEVEL in
%   absolute value.
%
%   STRICKCHECK differs from CLASSICCHECK() in that the tolerance on EPSLEVEL is
%   not relaxed by the length of the representation of F or by any finite
%   difference approximation of the gradient of F.
%
% See also STRICTCHECK, LOOSECHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with special cases ------------------------------------------------------

% Determine n (the length of the input):
n = length(f);

% Assume we're not happy:
ishappy = false; 

% Grab some preferences:
if ( nargin == 1 )
    pref = f.techPref();
    epslevel = pref.eps;
elseif ( isnumeric(pref) )
    epslevel = pref;
else
    epslevel = pref.eps;
end

% Convert scalar epslevel/tolerance inputs into vectors.
if ( isscalar(epslevel) )
    epslevel = repmat(epslevel, size(f.vscale));
end

% Deal with the trivial case:
if ( n < 2 )  % (Can't be simpler than a constant.)
    cutoff = n;
    return
end

% Check the vertical scale:
if ( isempty(vscl) || isempty(vscl) )
    if ( nargin < 2 || isempty(values) )
        % Compute some values if none were given:
        values = f.coeffs2vals(f.coeffs);
    end
    vscl = max(abs(values), [],  1);
else
    vscl = max(vscl, max(abs(values), [],  1));
end

if ( max(vscl) == 0 )
    % This is the zero function, so we must be happy.
    ishappy = true;
    cutoff = 1;
    return
elseif ( any(isinf(vscl(:))) )
    % Inf located. No cutoff.
    cutoff = n;
    epslevel = inf*epslevel;
    return
end

% NaNs are not allowed.
if ( any(isnan(f.coeffs(:))) )
    error('CHEBFUN:CHEBTECH:strictCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% Check for convergence and chop location --------------------------------------
testLength = min(n, max(5, round((n-1)/8))); 
ac = bsxfun(@rdivide, abs(f.coeffs), vscl);
f.coeffs(bsxfun(@le, ac, epslevel)) = 0;
tail = f.coeffs(end-testLength+1:end,:);
if ( ~any(tail(:)) )
    cutoff = find(max(f.coeffs, [], 2) > 0, 1, 'last');
    ishappy = true;
else
    cutoff = n;
    % Estimate the epslevel:
    epslevel = mean(abs(tail));
end

end
