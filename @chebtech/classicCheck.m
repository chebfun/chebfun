function [ishappy, epslevel, cutoff] = classicCheck(f, pref)
%CLASSICCHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = CLASSICCHECK(F) returns an estimated
%   location, the CUTOFF, at which the CHEBTECH F could be truncated to maintain
%   an accuracy of EPSLEVEL relative to F.VSCALE and F.HSCALE. ISHAPPY is
%   TRUE if CUTOFF < MIN(LENGTH(F.VALUES),2) or F.VSCALE = 0, and FALSE
%   otherwise.
%
%   [ISHAPPY, EPSLEVEL, CUTOFF] = CLASSICCHECK(F, PREF) allows additional
%   preferences to be passed. In particular, one can adjust the target accuracy
%   with PREF.CHEBTECH.EPS.
%
%   CLASSICCHECK first queries HAPPINESSREQUIREMENTS to obtain TESTLENGTH and
%   EPSLEVEL (see documentation below). If |F.COEFFS(1:TESTLENGTH)|/VSCALE <
%   EPSLEVEL, then the representation defined by F.VALUES and F.COEFFS is
%   deemed happy. The value returned in CUTOFF is essentially that from
%   TESTLENGTH (although it can be reduced if there are further COEFFS which
%   fall below EPSLEVEL).
%
%   HAPPINESSREQUIREMENTS defines what it means for a CHEBTECH to be happy.
%   [TESTLENGTH, EPSLEVEL] = HAPPINESSREQUIREMENTS(VALUES, COEFFS, VSCALE, PREF)
%   returns two scalars TESTLENGTH and EPSLEVEL. A CHEBTECH is deemed to be
%   'happy' if the coefficients COEFFS(1:TESTLENGTH) (recall that COEFFS are
%   stored in descending order) are all below EPSLEVEL. The default choice of
%   the test length is:
%       TESTLENGTH = n,             for n = 1:4
%       TESTLENGTH = 5,             for n = 5:44
%       TESTLENGTH = round((n-1)/8) for n > 44
%
%   EPSLEVEL is essentially the maximum of:
%       * eps*TESTLENGTH^(2/3)
%       * eps*grad (where grad is a finite difference approximation to the
%                   gradient of the function from VALUES. This is normalised by
%                   f.vscale and f.hscale).
%       * pref.chebtech.eps
%
%   Note that the accuracy check implemented in this function is the same as
%   that employed in Chebfun v4.x.
%
% See also STRICTCHECK, LOOSECHECK.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with special cases ------------------------------------------------------

% Determine n (the length of the input).
n = length(f);

% Assume we're not happy. (N'aww! :( )
ishappy = false;

% Grab some preferences:
if ( nargin == 1 )
    pref = f.pref();
    epslevel = pref.chebtech.eps;
elseif ( isnumeric(f) )
    epslevel = pref;
else
    epslevel = pref.chebtech.eps;
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

% NaNs are not allowed.
if ( any(isnan(f.coeffs(:))) )
    error('CHEBFUN:FUN:classicCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

% Check for convergence and chop location --------------------------------------

% Absolute value of coefficients, relative to vscale: (max across columns)
ac = max(bsxfun(@rdivide, abs(f.coeffs), f.vscale), [], 2);

% Take the maximum of the vscales:
vscale = max(f.vscale);

% Happiness requirements:
[testLength, epslevel] = ...
    happinessRequirements(f.values, f.coeffs, f.points(), vscale, f.hscale, epslevel);

if ( max(ac(1:testLength)) < epslevel )    % We have converged! Now chop tail:

    % We must be happy.
    ishappy = true;

    % Find first entry above epslevel:
    Tloc = find(ac >= epslevel, 1, 'first') - 1;

    % Check for the zero function!
    if ( isempty(Tloc) )
        cutoff = 1;
        return
    end

    % Compute the cumulative max of eps/4 and the tail entries:
    t = .25*eps;
    ac = ac(1:Tloc);               % Restrict to coefficients of interest.
    for k = 1:length(ac)           % Cumulative max.
        if ( ac(k) < t )
            ac(k) = t;
        else
            t = ac(k);
        end
    end

    % [TODO]: What does this mean?
    % Tbpb = Bang/buck of chopping at each pos:
    Tbpb = log(1e3*epslevel./ac) ./ (size(f.coeffs, 1) - (1:Tloc)');
    [ignored, Tchop] = max(Tbpb(3:Tloc));  % Tchop = position at which to chop.

    % We want to keep [c(0), c(1),  ..., c(cutoff)]:
    cutoff = n - Tchop - 2;

else

    % We're unhappy. :(
    cutoff = n;

end

end

function [testLength, epslevel] = ...
    happinessRequirements(values, coeffs, x, vscale, hscale, epslevel) %#ok<INUSL>
%HAPPINESSREQUIREMENTS   Define what it means for a CHEBTECH to be happy.
%   See documentation above.

% Grab the size:
n = size(values, 1);

% Length of tail to test.
testLength = min(n, max(5, round((n-1)/8)));

minPrec = 1e-4; % Worst case precision!

% Look at length of tail to loosen tolerance:
tailErr = min(minPrec, eps*testLength^(2/3));

% Look at finite difference gradient to loosen tolerance:
dx = diff(x)*ones(1, size(values, 2));
grad = (hscale/vscale) * norm(diff(values)./dx, inf);
gradErr = min(minPrec, eps*grad);

% Choose maximum between prescribed tolerance and estimated rounding errors:
epslevel = max([epslevel, gradErr, tailErr]);

end
