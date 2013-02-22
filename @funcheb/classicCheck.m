function [cutoff, epslevel] = classicCheck(f, pref)
%CLASSICCHECK  Attempt to trim trailing Chebyshev coefficients in a FUNCHEB.
%   [CUTOFF, EPSLEVEL] = CLASSICCHECK(F) returns an estimated location the
%   CUTOFF at which the FUNCHEB F could be truncated to maintain an accuracy of
%   EPSLEVEL relative to F.vscale and F.hscale.
%
%   [CUTOFF, EPSLEVEL] = CLASSICCHECK(F, PREF) allows additional preferences to
%   be passed. In particular, one can adjust the target accuracy with
%   PREF.FUNCHEB.EPS.
%
%   CLASSICCHECK first queries HAPPINESSREQUIREMENTS to obtain TESTLENGTH and
%   EPSLEVEL (see documentation below). If |F.COEFFS(1:TESTLENGTH)|/VSCALE <
%   EPSLEVEL, then the representation defined by F>VALUES and F.COEFFS is
%   deemed happy. The value returned in CUTOFF is essentially that from
%   TESTLENGTH (although it can be reduced if there are further COEFFS which
%   fall below EPSLEVEL).
%
%   HAPPINESSREQUIREMENTS defines what it means for a FUNCHEB to be happy.
%   [TESTLENGTH, EPSLEVEL] = HAPPINESSREQUIREMENTS(VALUES, COEFFS, VSCALE, PREF)
%   returns two scalars TESTLENGTH and EPSLEVEL. A FUNCHEB is deemed to be
%   'happy' if the coefficients COEFFS(1:TESTLENGTH) (recall that COEFFS are
%   stored in descending order) are all below EPSLEVEL. The default choice of
%   the test length is:
%       TESTLENGTH = n,             for n = 1:3
%       TESTLENGTH = 3,             for n = 4:25
%       TESTLENGTH = round((n-1)/8) for n > 25
%
%   EPSLEVEL is essentially the maximum of:
%       * eps*TESTLENGTH^(2/3)
%       * eps*grad (where grad is a finite difference approximation to the
%                   gradient of the function from VALUES. This is normalised by
%                   f.vscale and f.hscale).
%       * pref.funcheb.eps
%
% See also strictCheck.m, looseCheck.m.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with special cases -------------------------------------------------------

% Determine n (the length of the input).
n = length(f);

% Grab some preferences:
if ( nargin < 4 )
    pref = f.pref();
    epslevel = pref.(class(f)).eps;
elseif ( isnumeric(f) )
    epslevel = pref;
%     pref = f.pref();
end

% Deal with the trivial case 
if ( n < 2 ) % (Can't be simpler than a constant!)
    cutoff = n;
    return
end

% Check the vertical scale:
if ( max(f.vscale) == 0 )
    % This is the zero function, so we must be happy!
    cutoff = 1;
    return
elseif ( any(isinf(f.vscale)) )
    % Inf located. No cutoff.
    cutoff = n;
    return
end

% NaNs are not allowed
if ( any(isnan(f.coeffs(:))) )
    error('CHEBFUN:FUN:classicCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

% Check for convergence and chop location --------------------------------------

% Absolute value of coefficients, relative to vscale: (max across columns)
ac = max(bsxfun(@rdivide, abs(f.coeffs), f.vscale), [], 2);

% Take the minimum of the vscales:
vscale = max(f.vscale);

% Happiness requirements:
[testLength, epslevel] = ...
    happinessRequirements(f.values, f.coeffs, f.chebpts(length(f)), vscale, f.hscale, epslevel);

if ( max(ac(1:testLength)) < epslevel )    % We have converged! Now chop tail
    
    % Find first entry above epslevel:
    Tloc = find(ac >= epslevel, 1, 'first') - 1;
    
    % Check for the zero function!
    if ( isempty(Tloc) )                 
        cutoff = 1;
        return
    end
    
    % [TODO]: Figure out what the heck this does, and explain it.
    % Trim the coefficients:
%     ac = ac(1:Tloc);                      % Restrict to coeffs of interest
%     ac(1) = max(ac(1), .25*eps/vscale);   % Compute the cumulative max of 
%     for k = 2:Tloc                        %   eps/4/vscale and the tail entries
%         ac(k) = max(ac(k), ac(k-1));
%     end
    
    % Compute the cumulative max of eps/4 and the tail entries
    t = .25*eps;
    ac = ac(1:Tloc);                      % Restrict to coeffs of interest
    for k = 1:length(ac)                  % Cumulative max.
        if ( ac(k) < t )
            ac(k) = t;
        else
            t = ac(k);
        end
    end
    
    % Tbpb = Bang/buck of chopping at each pos:
    Tbpb = log(1e3*epslevel./ac) ./ (size(f.coeffs, 1) - (1:Tloc)');    
    [~, Tchop] = max(Tbpb(3:Tloc));       % Tchop = pos at which to chop

    % We want to keep [c(0), c(1),  ..., c(cutoff)]:
    cutoff = n - Tchop - 2;
    
else
    
    % We're unhappy. :(
    cutoff = n;
    
end

end

function [testLength, epslevel] = ...
    happinessRequirements(values, coeffs, x, vscale, hscale, epslevel) %#ok<INUSL>
%HAPPINESSREQUIREMENTS   Define what it means for a FUNCHEB to be happy.
%   See documentation above.

n = size(values, 1);
testLength = min(n, max(5, round((n-1)/8))); % Length of tail to test (see above)

minPrec = 1e-4; % Worst case precision! 

% Look at length of tail to loosen tolerance:
tailErr = min(minPrec, eps*testLength^(2/3));

% Look at finite difference gradient to loosen tolerance:
dx = max(diff(x), eps*hscale);
% dx = repmat(dx, 1, size(values, 2))
dx = dx*ones(1, size(values, 2));
grad = (hscale/vscale) * norm(diff(values)./dx, inf);
gradErr = min(minPrec, eps*grad);

% Choose maximum between prescribed tolerance and estimated rounding errors:
epslevel = max([epslevel, gradErr, tailErr]);

end





