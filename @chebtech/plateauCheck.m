function [ishappy, epsLevel, cutoff] = plateauCheck(f, values, pref)
%PLATEAUCHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = PLATEAUCHECK(F, VALUES) returns an estimated
%   location, the CUTOFF, at which the CHEBTECH F could be truncated. One of two
%   criteria must be met: Either:
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
%   [ISHAPPY, EPSLEVEL, CUTOFF] = PLATEAUCHECK(F, VALUES, PREF) allows additional
%   preferences to be passed. In particular, one can adjust the target accuracy
%   with PREF.EPS.
%
% See also STRICTCHECK, CLASSICCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab some preferences:
if ( nargin == 1 )
    pref = f.techPref();
    epslevel = pref.eps;
elseif ( isnumeric(pref) )
    epslevel = pref;
    pref = f.techPref();
else
    epslevel = pref.eps;
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
vscale = max( vscale, f.vscale );
absCoeff = absCoeff * diag(1./vscale);


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
% goal epslevel or the convergence appears to have levelled off for good
% (plateau).

% Guilty until proven innocent.
ishappy = false;
n = length(absCoeff);

%% 1. Strict test.

% Find the last place where the coeffs exceed the allowable level.
% Then go out a bit further to be safe.
cutoff = 4 + find( absCoeff >= epslevel, 1, 'last' );

if ( cutoff < 0.95*n )
    % Achieved the strict test.
    ishappy = true;
    
elseif ( n < 17 )
    % If there aren't enough coefficients, give up checking.
    cutoff = n;
    
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
    
    %%% Alternative windowed max: This avoids the for loop but might hog memory.
    %%index = bsxfun(@plus, (1:n)', 0:winsize-1);
    %%logabs = max(logabs(index),[],2);
    
    % Start with a low pass smoothing filter that introduces a lag.
    lag = 6;
    LPA = [1, zeros(1,lag-1), -2, zeros(1, lag-1), 1] / (lag^2);
    LPB = [1, -2, 1];
    smoothLAC = filter(LPA, LPB, winMax);  % smoothed logabs coeffs
    
    % If too little accuracy has been achieved, do nothing.
    tOK = find(smoothLAC < thresh, 1) - lag;
    if ( isempty(tOK) || (n - tOK < 16) )
        return
    end
    
    % Fit a least-squares line to windowed data.
    LSWindow = 24;
    A = [ ones(LSWindow,1), (1:LSWindow)'/LSWindow ];  % [1,x] quasimatrix
    % Create coefficients of all the windows: [ (1:LSW)', (2:LSW+1)', ... ]
    index = bsxfun(@plus, (0:LSWindow-1)', 1:length(smoothLAC)-LSWindow);
    LSLines = A \ smoothLAC(index);  % second row has all of the slopes
    slopes = LSLines(2,:);
        
    % Where is the decrease most rapid?
    slopeMin = min(slopes);
    
    % Don't look at anything until all substantial decrease has ended.
    tstart = find( slopes < 0.25*slopeMin, 1, 'last' );
    
    % Find where the decrease has permanently slowed to 1% of the fastest.
    isSlow = slopes(tstart:end) > 0.01*slopeMin;
    slow = tstart - 1 + find( isSlow );
    
    % Find the first run of 8 consecutive slow hits.
    first = find( slow(8:end) - slow(1:end-7) == 7, 1 );  % may be empty
    cutoff = slow(first);  % may be empty
    if ( isempty(cutoff) )
        cutoff = n;
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
