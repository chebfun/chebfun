function [ishappy, epslevel, cutoff] = plateauCheck(f, values, pref)
%PLATEAUCHECK   Attempt to trim trailing TRIGIER coefficients in a 
%TRIGTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = PLATEAUCHECK(F, VALUES) returns an 
%   estimated location, the CUTOFF, at which the TRIGTECH F could be 
%   truncated. One of two criteria must be met: Either:
%
%   (1) The coefficients are sufficiently small (as specified by the 
%   default EPS property of TRIGTECH) relative to F.VSCALE (or using 
%   absolute size if F.VSCALE=0); or
%
%   (2) The coefficients are somewhat small and apparently unlikely to
%   continue decreasing in a meaningful amount (i.e., have reached a 
%   "plateau" in convergence).
%
%   The reason for criterion (2) is that the problem may have a large 
%   condition number that prevents convergence to the full requested 
%   accuracy, as often happens in the collocation of differential 
%   equations.
%
%   Output EPSLEVEL is an estimate of the relative size of the last
%   "meaningful" expansion coefficients of the function, and the output 
%   CUTOFF is an estimate of how many of the coefficients are useful.
%
%   [ISHAPPY, EPSLEVEL, CUTOFF] = PLATEAUCHECK(F, VALUES, PREF) allows
%   additional preferences to be passed. In particular, one can adjust the
%   target accuracy with PREF.EPS.
%
% See also LINOPV4CHECK, STRICTCHECK, CLASSICCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab some preferences.
if ( nargin == 1 )
    pref = f.techPref();
    epslevel = pref.eps;
elseif ( isnumeric(pref) )
    epslevel = pref;
else
    epslevel = pref.eps;
end

% Grab the coefficients.
coeffs = f.coeffs;
N = length(coeffs);

% NaNs are not allowed.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:TRIGTECH:plateauCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% Set the function scaling for each vector of values.
maxvals = max(abs(values), [], 1);

% Check the vertical scale.
if ( max(maxvals) == 0 )
    % This is the zero function, so we must be happy!
    ishappy = true;
    cutoff = 1;
    return
elseif ( any(isinf(maxvals)) )
    % Inf located. No cutoff.
    ishappy = false;
    cutoff = N;
    return
end

% We add the (absolute value) of the coefficients correponsding to the same 
% degree together to use the same algorithm as the one used for CHEBTECH.
% Need to handle odd/even cases separately.
if ( mod(N, 2) == 0 ) % N even
    % [ -N/2 term; "positive" + "negative"; constant term ]
    absCoeffs = [ abs(coeffs(1,:)); ...
                  .5*abs(coeffs(2:N/2,:)) + .5*abs(coeffs(end:-1:N/2+2,:)); ...
                  abs(coeffs(N/2+1,:)) ];
else % N odd
    % [ "positive" + "negative"; constant term ]
    absCoeffs = ...
        [ .5*abs(coeffs(1:(N-1)/2,:)) + .5*abs(coeffs(end:-1:(N+3)/2,:)); ...
          abs(coeffs((N+1)/2,:)) ];
end
% Now N is odd in both cases.
N = size(absCoeffs, 1);

% We want a vector of coefficients with decreasing magnitude; again, to use
% the same algorithm as the one used for CHEBTECH.
absCoeffs = absCoeffs(end:-1:1,:);

% We omit the last 10% because aliasing can pollute them significantly.
N90 = ceil(0.90*N);
absCoeffs = absCoeffs(1:N90,:);
vscale = max(absCoeffs, [], 1);       
vscale = max([ vscale(:); f.vscale ]);
absCoeffs = absCoeffs/vscale;

numCol = size(absCoeffs, 2);
ishappy = false(1, numCol);
epslevels = zeros(1, numCol);
cutoff = zeros(1, numCol);
for m = 1:numCol
    [ishappy(m), epslevels(m), cutoff(m)] = ...
        checkColumn(absCoeffs(:,m), epslevel);
    if ( ~ishappy(m) )
        % No need to continue if it fails on any column.
        break
    end
end

epslevel = max(epslevels);
ishappy = all(ishappy); 

end

function [ishappy, epslevel, cutoff] = checkColumn(absCoeffs, epslevel)
% There are two ways to pass the test. Either the coefficients have 
% achieved the goal epslevel or the convergence appears to have levelled 
% off for good (plateau).

% Guilty until proven innocent.
ishappy = false;
N = length(absCoeffs);

%% 1. Strict test.

% Find the last place where the coeffs exceed the allowable level.
% Then go out a bit further to be safe.
cutoff = find(absCoeffs >= epslevel/50, 1, 'last');

if ( cutoff < 0.95*N )
    % Achieved the strict test.
    ishappy = true;
    
%% 2. Plateau test.
else
    
    % Demand at leastk=0 this much accuracy.
    thresh = max((2/3)*log(epslevel), log(1e-7));
    
    % Convergence is usually not far from linear in the log scale.
    logAbsCoeffs = log(absCoeffs);
    
    % Even though some methods can compute really small coefficients
    % relative to the norm, they ultimately contribute nothing. 
    % Also the occasional "zero" coefficient causes troublesome infinities.
    winMax = max(logAbsCoeffs, log(eps/1000));
    
    % Look for a sustained leveling off in the decrease.
    
    % [TODO]: Use the van Herk filter to do this more efficiently.
    
    % Start with a low pass smoothing filter that introduces a lag.
    lag = 6;
    LPA = [ 1, zeros(1,lag-1), -2, zeros(1, lag-1), 1 ]/(lag^2);
    LPB = [ 1, -2, 1 ];
    smoothLAC = filter(LPA, LPB, winMax);  % smoothed logabs coeffs
    
    % If too little accuracy has been achieved, do nothing.
    tOK = find(smoothLAC < thresh, 1) - lag;
    if ( isempty(tOK) || (N - tOK < 16) )
        return
    end
    
    % Fit a least-squares line to windowed data.
    LSWindow = 24;
    A = [ ones(LSWindow,1), (1:LSWindow)'/LSWindow ];  % [1,x] quasimatrix
    % Create coefficients of all the windows: [ (1:LSW)', (2:LSW+1)', ... ]
    index = bsxfun(@plus, (0:LSWindow-1)', 1:length(smoothLAC)-LSWindow);
    LSLines = A \ smoothLAC(index);  % second row has all of the slopes
    slopes = LSLines(2,:);
    slopes = [ nan(1,LSWindow), slopes ];
        
    % Where is the decrease most rapid?
    slopeMin = min(slopes);
    
    % Don't look at anything until all substantial decrease has ended.
    tstart = find(slopes < 0.3*slopeMin, 1, 'last');
    
    % Find where the decrease has slowed to 1% of the fastest.
    isSlow = slopes(tstart:end) > 0.01*slopeMin;
    slow = tstart - 1 + find(isSlow);
    
    % Find the first run of consecutive slow hits.
    numHits = 6;
    first = find(slow(numHits:end) - slow(1:end-numHits+1) == ... 
        numHits-1, 1);  % may be empty
    cutoff = slow(first);  % may be empty
    if ( isempty(cutoff) )
        cutoff = N;
    end
    
    % If the cut location is within the coefficient sequence, we're done.
    if ( cutoff < N )
        ishappy = true;
    end
    
end

% Deduce an epslevel. 
if ( ishappy == 1 )
    winEnd = min(N, cutoff + 6);
    epslevel = max(absCoeffs(cutoff:winEnd));
else
    epslevel = absCoeffs(N);
end
    
% Epslevel can't be better than eps:
epslevel = max(epslevel, eps);

end
