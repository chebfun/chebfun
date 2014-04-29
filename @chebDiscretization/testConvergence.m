function [isDone, epsLevel] = testConvergence(disc, values)
%TESTCONVERGENCE Happiness check.
%   Given a discretization, and a cell array of discretized functions,
%   check the equivalent Chebyshev polynomial representation for sufficient
%   convergence.

% We will test on an arbitrary linear combination of the individual
% functions.
s = 1 ./ (3*(1:numel(values))).';
newvalues = cell2mat(values(:).')*s;

% Convert to a piecewise chebfun.
u = toFunction(disc, newvalues);

% Test convergence on each piece. Start by obtaining the Chebyshev coefficients
% of all pieces, which we can then pass down to the testPiece method
coeffs = get(u, 'coeffs');
values = get(u, 'values');
d = disc.domain;
numInt = numel(d) - 1;
isDone = false(1, numInt);
epsLevel = 0;
for i = 1:numInt
    %    f = chebtech.constructor(values{i},u.vscale,hscale(i));
    [isDone(i),neweps] = plateauCheck(coeffs{i}, u.vscale);
    epsLevel = max( epsLevel, neweps );
end

end


function [ishappy, epslevel, cutoff] = plateauCheck(coeff, vscale)
%TODO: A summary documenting what's going on in this method would
% be nice. Also a short description of the outputs.

n = length(coeff);
epslevel = eps;

% Magnitude and rescale.
absCoeff = abs( coeff(end:-1:1) ) / vscale;

% NaNs are not allowed.
if ( any(isnan(coeff)) )
    error('CHEBFUN:FUN:plateauCheck:NaNeval', ...
        'Function returned NaN when evaluated.')
end

if ( vscale == 0 )
    % Trivially, the function is zero.
    ishappy = true;
    cutoff = 1;
    return
end


%% Serious checking starts here.

% The strict test is for the coefficients to get and stay below epslevel.

% Starting from the tail end, where do the coefficients first exceed the
% allowable level?
cutoff = find( absCoeff(end:-1:1) > epslevel, 1 );

% How many do we require in the tail?
testLength = min(n, max(5, round((n-1)/8))); 


if ( cutoff > testLength)
    %% Strict check passed.
    ishappy = true;
    
elseif ( n < 17 )
    %% If there aren't enough coefficients, give up checking.
    ishappy = false;
    epslevel = absCoeff(n);
    cutoff = n;
    
    
else
    %% Apply the plateau test.
    ishappy = false;
    thresh = 1e-6;  % demand at least this much accuracy
    
    % Convergence is usually not far from linear in the log scale.
    logabs = log(absCoeff);
    
    % Even though ultraS can compute really small coefficients relative to the
    % norm, they ultimately contribute nothing. Also the occasional "zero"
    % coefficient causes troublesome infinities.
    logabs = max( logabs, log(eps/10) );
    
    % Goal: Look for a sustained leveling off in the decrease of the
    % coefficients.
    
    % Start with a low pass filter that introduces a lag.
    lag = 8;
    LPA = [1 zeros(1,lag-1) -2 zeros(1,lag-1) 1]/(lag^2);
    LPB = [1 -2 1];
    smoothLAC = filter( LPA, LPB, logabs );
    
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
    tstart = find( smoothDiff < 0.4*SDmin, 1, 'last' );
    
    % Find where the decrease has slowed to 10% of fastest.
    slow = tstart - 1 + find( smoothDiff(tstart:end) > 0.1*SDmin );
    slow = slow - floor(lag/2);  % partly compensate for the filter lag
    
    % Find the first run of 5 consecutive slow hits.
    first = find( slow(5:end) - slow(1:end-4) == 4, 1 );  % may be empty
    cutoff = slow(first);  % may be empty
    
    % If the cut location is within the coefficient sequence, we're done.
    if ( cutoff < n )
        ishappy = true;
        % Use the information from the cut to deduce an eps level.
        window = min( n, cutoff+(1:4) );
        epslevel = exp( max( logabs(window) ) );
    end
    
end

end

