function [isDone, epsLevel] = testConvergence(disc, values, vscale, pref)
%TESTCONVERGENCE Happiness check.
%   Given a discretization, and a cell array of discretized functions,
%   check the equivalent Chebyshev polynomial representation for sufficient
%   convergence.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

if ( nargin < 4 )
    pref = cheboppref;
    if ( nargin < 3 )
        vscale = Inf;   % will have no effect
    end
end

% We will test on an arbitrary linear combination of the individual
% functions.
s = 1 ./ (3*(1:numel(values))).';
newValues = cell2mat(values(:).')*s;

% Convert to a piecewise chebfun.
u = toFunctionOut(disc, newValues);

% Test convergence on each piece. Start by obtaining the Chebyshev coefficients
% of all pieces, which we can then pass down to the testPiece method
coeffs = get(u, 'coeffs');
d = disc.domain;
numInt = numel(d) - 1;
isDone = false(1, numInt);
epsLevel = 0;

% If an external vscale was supplied, it can supplant the inherent scale of
% the result.
vscale = min(u.vscale,vscale);

for i = 1:numInt
    [isDone(i),neweps] = plateauCheck(coeffs{i}, vscale, pref);
    epsLevel = max( epsLevel, neweps );
end

end

function [ishappy, epslevel, cutoff] = plateauCheck(coeff, vscale, pref)
%  Inputs: 
%    coeff: vector of Cheb poly coefficients (high order to low)
%    vscale: indication of the scale to resolve relative to (default=Inf,
%           no effect)
%    pref: cheboppref
%
%  Outputs:
%    ishappy: true if convergence was achieved
%    epslevel: the apparent epslevel of the truncation
%    cutoff: where to truncate the coefficients

% This check is needed because of condition numbers in differential
% equations. We can't be sure that a solution will ever be resolved to full
% precision, so we have to be willing to stop if the convergence appears to
% have trailed off.

% NaNs are not allowed.
if ( any(isnan(coeff)) )
    error('CHEBFUN:FUN:plateauCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

n = length(coeff);

% We omit the last 12% because aliasing can pollute them significantly.
n = ceil( 0.88*n );

epslevel = pref.errTol;  % preferred tolerance

% Magnitude and rescale.
absCoeff = abs( coeff(n:-1:1) ) / vscale;

if ( vscale == 0 )
    % Trivially, the function is zero.
    ishappy = true;
    cutoff = 1;
    return
end

%% Serious checking starts here.
% There are two ways to pass the test. Either the coefficients have
% achieved the goal epslevel, or the convergence appears to have levelled
% off for good (plateau). 

%TODO: Unify and locate with the chebtech happiness checks. 

% 1. Strict test.

% Find the last place where the coeffs exceed the allowable level.
% Then go out a bit further to be safe.
cutoff = 4 + find( absCoeff >= epslevel, 1, 'last' );

if ( cutoff < 0.95*n ) 
    % Achieved the strict test. 
    ishappy = true;
    
elseif ( n < 17 )
    % If there aren't enough coefficients, give up checking.
    ishappy = false;
    epslevel = absCoeff(n);
    cutoff = n;
    
% 2. Plateau test.
else
    ishappy = false;
    
    % Demand at least this much accuracy.
    thresh = max(log(epslevel),log(1e-7));  
    
    % Convergence is usually not far from linear in the log scale.
    logabs = log(absCoeff);
    
    % Even though some methods can compute really small coefficients relative to the
    % norm, they ultimately contribute nothing. Also the occasional "zero"
    % coefficient causes troublesome infinities.
    logabs = max( logabs, log(eps/1000) );
    
    % Look for a sustained leveling off in the decrease.
    
    % Symmetries can cause one or more consecutive coefficients to be zero,
    % and we only care about the nonzero ones. Use a windowed max to remove
    % the small values.
    % TODO: Use the van Herk filter to do this more efficiently.
    winsize = 6;
    winmax = logabs;
    for k = 1:winsize
        winmax = max(winmax(1:end-1),logabs(k+1:end));
    end
    n = length(winmax);

    %%% Alternative windowed max: This avoids the for loop but might hog memory.
    %%index = bsxfun(@plus, (1:n)', 0:winsize-1);
    %%logabs = max(logabs(index),[],2);
    
    % Start with a low pass smoothing filter that introduces a lag.
    lag = 6;
    LPA = [1, zeros(1,lag-1), -2, zeros(1, lag-1) 1]/(lag^2);
    LPB = [1, -2, 1];
    smoothLAC = filter( LPA, LPB, winmax );  % smoothed logabs coeffs
    
    % If too little accuracy has been achieved, do nothing.
    tOK = find(smoothLAC < thresh, 1) - lag;
    if ( isempty(tOK) || n - tOK < 16 )
        return
    end
    
    % Smooth the first difference of the smoothed coefficent sequence.
    smoothDiff = filter( LPA, LPB, diff(smoothLAC) );
    
    % Where is the decrease most rapid?
    SDmin = min(smoothDiff);
    
    % Don't look at anything until all substantial decrease has ended.
    tstart = find( smoothDiff < 0.25*SDmin, 1, 'last' );
    
    % Find where the decrease has permanently slowed to 10% of the fastest.
    isSlow = smoothDiff(tstart:end) > 0.01*SDmin;
    lastFast = find(~isSlow,1,'last');
    if ( isempty(lastFast) )
        lastFast = 0;
    end
    slow = tstart + lastFast - 1 + find( isSlow(lastFast+1:end) );
    slow = slow - floor(lag/2);  % compensate for the filter lag
    
    % Find the first run of 5 consecutive slow hits.
    first = find( slow(5:end) - slow(1:end-4) == 4, 1 );  % may be empty
    cutoff = slow(first);  % may be empty, will give false next
    
    % If the cut location is within the coefficient sequence, we're done.
    if ( cutoff < n )
        ishappy = true;
    end
    
end

if ( ishappy )
    % Use the information from the cut to deduce an eps level.
    winEnd = min( n, cutoff + 4 );
    epslevel = max( absCoeff(cutoff:winEnd) );
end

end
