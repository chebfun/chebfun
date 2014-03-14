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
numInt = numel(disc.domain) - 1;
isDone = true(1, numInt);
epsLevel = 0;
for i = 1:numInt
    [isDone(i), t2] = testPiece(coeffs{i}.');
    epsLevel = max(epsLevel, t2);
end

end


function [isDone, epsLevel] = testPieceOld(coeffs)
% Test convergence on a single subinterval.

% FIXME: This is the v4 test. It still has an ad hoc nature.

isDone = false;
epsLevel = eps;
thresh = 1e-6;  % demand at least this much accuracy

% Flip the coeffs to get zero degree first
c = coeffs(end:-1:1);

n = length(c);
if n < 30, return, end

% Magnitude and rescale.
ac = abs(c);
ac = ac/min(max(ac), 1);

% Smooth using a windowed max to dampen symmetry oscillations.
maxac = ac;
for k = 1:8
    maxac = max(maxac(1:end - 1), ac(k + 1:end));
end

% If too little accuracy has been achieved, do nothing.
t = find(maxac < thresh, 1);
if ( isempty(t) || n - t < 16 )
    return
end

% Find where improvement in the windowed max seems to stop, by looking at
% the derivative of a smoother form of the curve.
dmax = diff( conv( [1 1 1 1]/4, log(maxac(t:end)) ) );
mindmax = dmax;
for k = 1:2
    mindmax = min(mindmax(1:end - 1), dmax(k + 1:end));
end

% Find the place to chop the series.
cutHere = find(mindmax > 0.01*min(mindmax), 3);
if ( isempty(cutHere) )
    cutHere = 1;
else
    cutHere = cutHere(end) + t + k + 3;
end

% If the cut location is within the given coefficients, we're done.
if ( cutHere < n )
    isDone = true;
    % Use the information from cut to deduce an epsLevel
    epsLevel = max( abs(c(cutHere + 1)) );
end

end

function [isDone, epsLevel] = testPiece(coeffs)
% Test convergence on a single subinterval.

isDone = false;
epsLevel = eps;
thresh = 1e-6;  % demand at least this much accuracy

% Flip the coeffs to get zero degree first
c = coeffs(end:-1:1);

n = length(c);
if n < 30, return, end

% Magnitude and rescale.
absC = abs(c);
absC = absC/min(max(absC), 1);

% If full accuracy has been achieved, we will stop. Check for 10
% consecutive coefficients well under machine eps.
% This can happen with ultraS, for example.
tiny = find( absC < eps/1000 );
if any( tiny(10:end) - tiny(1:end-9) == 9 )
    isDone = true;
    epsLevel = eps;
    
else
    % Look for a sustained leveling off in the decrease of the
    % coefficients.
    logabsC = log(absC);
    logabsC = max( logabsC, log(eps/10) );
    
    % Use a low pass filter.
    lag = 8;
    LPA = [1 zeros(1,lag-1) -2 zeros(1,lag-1) 1]/(lag^2);
    LPB = [1 -2 1];
    smoothLAC = filter( LPA, LPB, logabsC );
    
    % If too little accuracy has been achieved, do nothing.
    tOK = find(smoothLAC < thresh, 1) - lag;
    if ( isempty(tOK) || n - tOK < 16 )
        return
    end
    
    
    % To use rate of decrease, smooth the first difference of the sequence.
    smoothDiff = filter( LPA, LPB, diff(smoothLAC) );
    
    % Where is the decrease most rapid?
    [SDmin,tmin] = min(smoothDiff);
    tmin = tmin - lag + 1;
    
    % Find where decrease has slowed dramatically.
    slow = tmin - 1 + find( smoothDiff(tmin:end) > 0.05*SDmin );
    
    % Find the first place with 5 consecutive slow hits.
    first = find( slow(5:end) - slow(1:end-4) == 4, 1 );  % may be empty
    cutHere = slow(first);  % may be empty
    
    % If the cut location is within the given coefficients, we're done.
    if ( cutHere < n )
        isDone = true;
        % Use the information from the cut to deduce an epsLevel.
        window = min( n, cutHere+(1:4) );
        epsLevel = exp( max( logabsC(window) ) );
    end
end
end
