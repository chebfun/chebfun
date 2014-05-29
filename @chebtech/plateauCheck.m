function [ishappy, epslevel, cutoff] = plateauCheck(f, values, pref)
%PLATEAUCHECK   Seek a plateau in Chebyshev coefficients.
%  Inputs:
%    f:      chebtech
%    vscale: indication of the scale to resolve relative to (default=Inf,
%            no effect)
%    pref:   cheboppref
%
%  Outputs:
%    ishappy:  true if convergence was achieved
%    epslevel: the apparent epslevel of the truncation
%    cutoff:   where to truncate the coefficients
%
% This check is needed because of condition numbers in differential equations.
% We can't be sure that a solution will ever be resolved to full precision, so
% we have to be willing to stop if the convergence appears to have trailed off.

coeff = f.coeffs;
vscale = 0;
if ( nargin > 1 ) 
    if isnumeric(values)
        vscale = norm(values(:),inf);
    else
        pref = values;
    end
else
    pref = f.techpref();
end

% NaNs are not allowed.
if ( any(isnan(coeff)) )
    error('CHEBFUN:FUN:plateauCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

% We omit the last 12% because aliasing can pollute them significantly.
n = length(coeff);
n88 = ceil( 0.88*n );
% Preferred tolerance
epslevel = pref.eps;  
% Magnitude and rescale.
if ( vscale > 0 )
    absCoeff = abs( coeff(n:-1:1) ) / vscale;
end

% %%%%%%%%%%%%%%%%%%%%%%%% Serious checking starts here. %%%%%%%%%%%%%%%%%%%%%%%
% There are two ways to pass the test. Either the coefficients have
% achieved the goal epslevel, or the convergence appears to have levelled
% off for good (plateau).

% Guilty until proven innocent.
ishappy = false;

%% 1. Strict test.

% Find the last place where the coeffs exceed the allowable level.
% Then go out a bit further to be safe.
cutoff = 4 + find( absCoeff >= epslevel, 1, 'last' );

if ( cutoff < 0.95*n88 )
    % Achieved the strict test.
    ishappy = true;
    
elseif ( n88 < 17 )
    % If there aren't enough coefficients, give up checking.
    epslevel = absCoeff(n88);
    cutoff = n88;
    
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
    
    % TODO: Use the van Herk filter to do this more efficiently.
    
    % Symmetries can cause one or more consecutive coefficients to be zero, and
    % we only care about the nonzero ones. Use a windowed max to remove the
    % small values.
    winSize = 6;
    winMax = logAbs;
    for k = 1:winSize
        winMax = max( winMax(1:end-1), logAbs(k+1:end) );
    end
    n88 = length(winMax);
    
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
    if ( isempty(tOK) || (n88 - tOK < 16) )
        return
    end
    
    % Smooth the first difference of the smoothed coefficent sequence.
    smoothDiff = filter( LPA, LPB, diff(smoothLAC) );
    
    % Where is the decrease most rapid?
    SDmin = min(smoothDiff);
    
    % Don't look at anything until all substantial decrease has ended.
    tstart = find( smoothDiff < 0.25*SDmin, 1, 'last' );
    
    % Find where the decrease has permanently slowed to 1% of the fastest.
    isSlow = smoothDiff(tstart:end) > 0.01*SDmin;
    lastFast = find(~isSlow, 1, 'last');
    if ( isempty(lastFast) )
        lastFast = 0;
    end
    slow = tstart + lastFast - 1 + find( isSlow(lastFast+1:end) );
    slow = slow - floor(lag/2);  % compensate for the filter lag
    
    % Find the first run of 5 consecutive slow hits.
    first = find( slow(5:end) - slow(1:end-4) == 4, 1 );  % may be empty
    cutoff = slow(first);  % may be empty, will give false next
    
    % If the cut location is within the coefficient sequence, we're done.
    if ( cutoff < n88 )
        ishappy = true;
    end
    
end

if ( ishappy )
    % Use the information from the cut to deduce an eps level.
    winEnd = min( n88, cutoff + 4 );
    epslevel = max( absCoeff(cutoff:winEnd) );
end

end
