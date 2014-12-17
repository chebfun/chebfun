function [ishappy, epslevel, cutoff] = classicCheck(f, pref)
%CLASSICCHECK   Attempt to trim trailing Fourier coefficients in a TRIGTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = CLASSICCHECK(F, VALUES) returns an estimated
%   location, the CUTOFF, at which the TRIGTECH F could be truncated to maintain
%   an accuracy of EPSLEVEL relative to F.VSCALE and F.HSCALE. ISHAPPY is TRUE
%   if CUTOFF < MIN(LENGTH(VALUES),2) or F.VSCALE = 0, and FALSE otherwise.
%   If ISHAPPY is false, EPSLEVEL returns an estimate of the accuracy achieved.
%
%   [ISHAPPY, EPSLEVEL, CUTOFF] = CLASSICCHECK(F, PREF) allows additional
%   preferences to be passed. In particular, one can adjust the target accuracy
%   with PREF.EPS.
%
%   CLASSICCHECK first queries HAPPINESSREQUIREMENTS to obtain TESTLENGTH and
%   EPSLEVEL (see documentation below). If |F.COEFFS(1:TESTLENGTH)|/VSCALE <
%   EPSLEVEL, then the representation defined by F.COEFFS is deemed happy. The
%   value returned in CUTOFF is essentially that from TESTLENGTH (although it
%   can be reduced if there are further COEFFS which fall below EPSLEVEL).
%
%   HAPPINESSREQUIREMENTS defines what it means for a TRIGTECH to be happy.
%   [TESTLENGTH, EPSLEVEL] = HAPPINESSREQUIREMENTS(VALUES, COEFFS, VSCALE, PREF)
%   returns two scalars TESTLENGTH and EPSLEVEL. A TRIGTECH is deemed to be
%   'happy' if the coefficients COEFFS(1:TESTLENGTH) (recall that COEFFS are
%   stored in ascending order) are all below EPSLEVEL. The default choice of
%   the test length is:
%       TESTLENGTH = n,             for n = 1:2
%       TESTLENGTH = 3,             for n = 3:44
%       TESTLENGTH = round((n-1)/8) for n > 44
%
%   EPSLEVEL is essentially the maximum of:
%       * pref.eps
%       * eps*TESTLENGTH
%       * eps*condEst (where condEst is an estimate of the condition number
%                      based upon a finite difference approximation to the
%                      gradient of the function from F.VALUES.).
%   However, the final two estimated values can be no larger than 1e-4.
%
%   Note that the accuracy check implemented in this function is the (roughly)
%   same as that employed in Chebfun v4.x.
%
% See also STRICTCHECK, LOOSECHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with special cases ------------------------------------------------------

% Determine n (the length of the input).
n = length(f);

% Assume we're not happy. (N'aww! :( )
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

% If one column of f is the zero function, we will get into trouble further
% down when we take the absolute value of the coefficients relative to vscale
% and compute the relative condition number estimate in happinessCheck.  We
% replace zero vscales by eps to avoid division by zero.
ind = f.vscale == 0;
f.vscale(ind) = eps;

% NaNs are not allowed.
if ( any(isnan(f.coeffs(:))) )
    error('CHEBFUN:FUN:classicCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

% Do the test on the vector formed by the sum of the absolute value of the
% positive and negative mode coefficients.
% [TODO] This reversal of coeffs can be removed but then classicCheck will need
% to be written carefully:
absCoeffs = abs(f.coeffs(end:-1:1,:));

% Need to handle odd/even cases separately.
isEven = ~mod(n,2);
if isEven
    % In this case the negative cofficients have an additional term
    % corresponding to the cos(N/2*x) coefficient.
    f.coeffs = [absCoeffs(n,:);absCoeffs(n-1:-1:n/2+1,:)+absCoeffs(1:n/2-1,:);absCoeffs(n/2,:)];
else
    f.coeffs = [absCoeffs(n:-1:(n+1)/2+1,:)+absCoeffs(1:(n+1)/2-1,:);absCoeffs((n+1)/2,:)];
end

n = size(f.coeffs, 1);

% Check for convergence and chop location --------------------------------------

% Absolute value of coefficients, relative to vscale: (max across columns)
ac = bsxfun(@rdivide, abs(f.coeffs), f.vscale);

% Happiness requirements:
[testLength, epslevel] = ...
    happinessRequirements(f.values, f.coeffs, f.points(), f.vscale, f.hscale, epslevel);

if ( all(max(ac(1:testLength, :)) < epslevel) ) % We have converged! Chop tail:
    % We must be happy.
    ishappy = true;

    % Find first row of coeffs with entry above epslevel:
    rowsWithLargeCoeffs = any(bsxfun(@ge, ac, epslevel), 2);
    Tloc = find(rowsWithLargeCoeffs, 1, 'first') - 1;

    % Check for the zero function!
    if ( isempty(Tloc) )
        cutoff = 1;
        return
    end

    % Compute the cumulative max of eps/4 and the tail entries:
    t = .25*eps*ones(1, size(ac, 2));
    ac = ac(1:Tloc, :);             % Restrict to coefficients of interest.
    for k = 1:size(ac, 1)           % Cumulative maximum.
        ind = ac(k, :) < t;
        ac(k, ind) = t(ind);
        ind = ac(k, :) >= t;
        t(ind) = ac(k, ind);
    end

    % Obtain an estimate for much accuracy we'd gain compared to reducing
    % length ("bang for buck"):
    bang = log(1e3*bsxfun(@rdivide, epslevel, ac));
    buck = n - (1:Tloc).';
    Tbpb = bsxfun(@rdivide, bang, buck);

    % Compute position at which to chop.  Keep greatest number of coefficients
    % demanded by any of the columns.
    [ignored, perColTchop] = max(Tbpb(3:Tloc, :));
    Tchop = min(perColTchop);

    % Cutoff value.
    cutoff = n - Tchop - 2;
    % We want to keep [c(-cutoff), ...,c(-1), c(0), c(1), ..., c(cutoff)],
    % so the number of coefficients that will be thrown away is actually 2*cutoff-1.
    cutoff = 2*cutoff-1;
    
else

    % We're unhappy. :(
    cutoff = n;

end

end

function [testLength, epslevel] = ...
    happinessRequirements(values, coeffs, x, vscale, hscale, epslevel) %#ok<INUSL>
%HAPPINESSREQUIREMENTS   Define what it means for a TRIGTECH to be happy.
%   See documentation above.

% Grab the size:
n = size(coeffs, 1);

% We will not allow the estimated rounding errors to be cruder than this value:
minPrec = 1e-4; % Worst case precision!

% Length of tail to test.
testLength = min(n, max(3, round((n-1)/8)));

% Look at length of tail to loosen tolerance:
tailErr = eps*testLength;
tailErr = min(tailErr, minPrec);

% Estimate the condition number of the input function by
% ||f(x+eps(x)) - f(x)||_inf / ||f||_inf ~~ (eps(hscale)/vscale)*f'.
dy = diff(values);
dx = diff(x)*ones(1, size(values, 2));
gradEst = max(abs(dy./dx));             % Finite difference approx.
condEst = eps(hscale)./vscale.*gradEst; % Condition number estimate.
condEst = min(condEst, minPrec);      

% Choose maximum between prescribed tolerance and estimated rounding errors:
epslevel = max(max(epslevel, condEst), tailErr);

end
