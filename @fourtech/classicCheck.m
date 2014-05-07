function [ishappy, epslevel, cutoff] = classicCheck(f, pref)
%CLASSICCHECK   Attempt to trim trailing Fourier coefficients in a FOURTECH.
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
elseif ( isnumeric(f) )
    epslevel = pref;
else
    epslevel = pref.eps;
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

% Do the test on the vector formed by the sum of the absolute value of the
% positive and negative mode coefficients.
absCoeffs = abs(f.coeffs);

% Need to handle odd/even cases separately
isEven = ~mod(n,2);
if isEven
    % In this case the positive cofficients have an additional term
    % corresponding to the cos(N/2*x) coefficient.
%     f.coeffs = [absCoeffs(n/2,:);absCoeffs(n/2+1:n-1,:)+absCoeffs(n/2-1:-1:1,:);absCoeffs(n,:)];
    f.coeffs = [absCoeffs(n,:);absCoeffs(n-1:-1:n/2+1,:)+absCoeffs(1:n/2-1,:);absCoeffs(n/2,:)];
else
%     f.coeffs = [absCoeffs((n+1)/2,:);absCoeffs((n+1)/2+1:n,:)+absCoeffs((n+1)/2-1:-1:1,:)];
    f.coeffs = [absCoeffs(n:-1:(n+1)/2+1,:)+absCoeffs(1:(n+1)/2-1,:);absCoeffs((n+1)/2,:)];
end

% c = abs(f.coeffs);
% m = floor(n/2);
% c2 = c(1:ceil(n/2),:);
% c2(1:m,:) = c2(1:m,:) + flipud(c(end-m+1:end,:));
% f.coeffs = c2;
n = size(f.coeffs, 1);

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
    for k = 1:length(ac)           % Cumulative maximum.
        if ( ac(k) < t )
            ac(k) = t;
        else
            t = ac(k);
        end
    end

    % Obtain an estimate for much accuracy we'd gain compared to reducing
    % length ("bang for buck"):
    Tbpb = log(1e3*epslevel./ac) ./ (n - (1:Tloc)');
    [~, Tchop] = max(Tbpb(3:Tloc));  % Position at which to chop.

    % We want to keep [c(0), c(1), ..., c(cutoff)]:
    cutoff = n - Tchop - 2;

else

    % We're unhappy. :(
    cutoff = n;

end

end

function [testLength, epslevel] = ...
    happinessRequirements(values, coeffs, x, vscale, hscale, epslevel) %#ok<INUSL>
%HAPPINESSREQUIREMENTS   Define what it means for a FOURTECH to be happy.
%   See documentation above.

% Grab the size:
n = size(coeffs, 1);

% We will not allow the estimated rounding errors to be cruder than this value:
minPrec = 1e-4; % Worst case precision!

% Length of tail to test.
% testLength = min(n, max(5, round((n-1)/8)));
testLength = min(n, max(4, round((n-1)/8)));

% Look at length of tail to loosen tolerance:
tailErr = eps*testLength^(2/3);
tailErr = min(tailErr, minPrec);

% Estimate the condition number of the input function by
%    ||f(x+eps(x)) - f(x)||_inf / ||f||_inf ~~ (eps(hscale)/vscale)*f'.
dy = diff(values);
dx = diff(x)*ones(1, size(values, 2));
gradEst = norm(dy./dx, inf);          % Finite difference approx.
condEst = eps(hscale)/vscale*gradEst; % Condition number estimate.
condEst = min(condEst, minPrec);      

% Choose maximum between prescribed tolerance and estimated rounding errors:
epslevel = max([epslevel, tailErr, condEst]);

end
