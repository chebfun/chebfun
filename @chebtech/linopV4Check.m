function [ishappy, epsLevel, cutoff] = linopV4Check(f, values, pref)
%LINOPV4CHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = LINOPV4CHECK(F, VALUES) returns an estimated
%   location, the CUTOFF, at which the CHEBTECH F could be truncated. It's
%   based on the same functionality provided in Version 4 of Chebfun and is
%   more aggressive about truncation than the alternative PLATEAUCHECK.
%
%   The first output measures "happiness" (sufficient resolution). One of two
%   criteria must be met. Either:
%
%     (1) The coefficients are sufficiently small (as specified by the default
%     EPS property of CHEBTECH) relative to F.VSCALE (or using absolute size if
%     F.VSCALE=0); or
%
%     (2) The coefficients are somewhat small and apparently unlikely to
%     continue decreasing in a meaningful amount (i.e., have reached a "plateau"
%     in convergence).
%
%   The reason for criterion (2) is that the problem may have a large condition
%   number that prevents convergence to the full requested accuracy, as often
%   happens in the collocation of differential equations.
%
%   Output EPSLEVEL is an estimate of the relative size of the last
%   "meaningful" expansion coefficients of the function, and the output 
%   CUTOFF is an estimate of how many of the coefficients are useful.
%
%   [...] = LINOPV4CHECK(F, VALUES, PREF) allows additional
%   preferences to be passed. In particular, one can adjust the target 
%   accuracy with PREF.EPS.
%
% See also PLATEAUCHECK, STRICTCHECK, CLASSICCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab some preferences:
if ( nargin == 1 )
    pref = f.techPref();
    epsLevel = pref.eps;
elseif ( isnumeric(pref) )
    epsLevel = pref;
    pref = f.techPref();
else
    epsLevel = pref.eps;
end

% Grab the coefficients:
coeff = f.coeffs;
n = length(coeff);

% NaNs are not allowed.
if ( any(isnan(coeff)) )
    error('CHEBFUN:FUN:plateauCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

% Set the function scaling for each vector of values.
maxvals = max(abs(values), [], 1);

% Check the vertical scale:
if ( max(maxvals) == 0 )
    % This is the zero function, so we must be happy!
    ishappy = true;
    cutoff = 1;
    return
elseif ( any(isinf(maxvals)) )
    % Inf located. No cutoff.
    cutoff = n;
    return
end

%%
% We omit the last 10% because aliasing can pollute them significantly.
n90 = ceil( 0.90*n );
absCoeff = abs( coeff(end:-1:end+1-n90,:) );  % switch to low->high ordering
vscale = max(absCoeff,[],1);          % scaling in each column
vscale = max( [vscale(:); f.vscale] );
absCoeff = absCoeff / vscale;


%% Deal with array-valued functions.

numCol = size(coeff, 2);
ishappy = false(1,numCol);
epsLevel = zeros(1,numCol);
cutoff = zeros(1,numCol);
for m = 1:numCol
    [ishappy(m), epsLevel(m), cutoff(m)] = checkColumn(absCoeff(:,m),pref.eps);
    if ( ~ishappy(m) )
        % No need to continue if it fails on any column.
        break
    end
end

epsLevel = max(epsLevel);
ishappy = all(ishappy); 

end

%% %%%%%%%%%%%%%%%%%%%%%%% Serious checking starts here. %%%%%%%%%%%%%%%%%%%%%%%
function [ishappy, epslevel, cutoff] = checkColumn(absCoeff,epslevel)

% There are two ways to pass the test. Either the coefficients have achieved the
% goal epslevel or the convergence appears to have leveled off for good

% Guilty until proven innocent.
ishappy = false;
n = length(absCoeff);

%% 1. Strict test.

% Find the last place where the coeffs exceed the allowable level.
% Then go out a bit further to be safe.
cutoff = 4 + find( absCoeff >= epslevel/50, 1, 'last' );

if ( cutoff < 0.95*n )
    % Achieved the strict test.
    ishappy = true;
        
%% 2. Plateau test.
else
    
    % Demand at least this much accuracy.
    thresh = max((2/3)*log(epslevel), log(1e-7));
    
    
    
    % Convergence is usually not far from linear in the log scale.
    logAbs = log(absCoeff);
    
    % Even though some methods can compute really small coefficients relative to
    % the norm, they ultimately contribute nothing. Also the occasional "zero"
    % coefficient causes troublesome infinities.
    logAbs = max( logAbs, log(eps/1000) );
    
    % Look for a sustained leveling off in the decrease.
    
    % TODO: Use the van Herk filter to do this more efficiently. (Low priority.)
    
    % Symmetries can cause one or more consecutive coefficients to be zero, and
    % we only care about the nonzero ones. Use a windowed max to remove the
    % small values.
    winSize = 6;
    winMax = logAbs;
    for k = 1:winSize
        winMax = max( winMax(1:end-1), logAbs(k+1:end) );
    end
    n = length(winMax);
    
    % If too little accuracy has been achieved, do nothing.
    tOK = find(winMax < thresh, 1);
    if ( isempty(tOK) || (n - tOK < 16) )
        return
    end
    
    % Smooth, then difference the windowed max.
    diffMax = diff( conv( [1 1 1 1]/4, winMax(tOK:end) ) );

    % Do a windowed min. 
    winMinDiffMax = diffMax;
    winMinDiffMax = min(winMinDiffMax(1:end-1),diffMax(2:end));

    cutoff = find(winMinDiffMax > 0.01*min(winMinDiffMax), 3);
    if ( isempty(cutoff) )
        cutoff = n;
    else
        cutoff = cutoff(end) + tOK + 2 + 3;
    end
   
    % If the cut location is within the coefficient sequence, we're done.
    if ( cutoff < n )
        ishappy = true;
    end
    
end

% Deduce an epslevel. 
if ( ishappy )
    winEnd = min( n, cutoff + 4 );
    epslevel = max( absCoeff(cutoff:winEnd) );
else
    epslevel = absCoeff(n);
end
    
% Epslevel can't be better than eps:
epslevel = max(epslevel, eps);

end
