function poleOrder = findPoleOrder(op, singEnd)
%FINDPOLEORDER   Finds the order of the pole in the function handle OP at
%   x = 1 or -1 depending upon the string 'left' or 'right' passed in 
%   SINGEND.
%   
% Example:
%   p = singfun.findPoleOrder(@(x) 1./(1 - x), 'right')
%   p = singfun.findPoleOrder(@(x) 1./(1 + x).^2, 'left')
%
% See also FINDSINGORDER, FINDSINGEXPONENTS.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% The algorithm tries to find the smallest non-negative integer
% k such that op(x).*(1-x)^k is bounded at 1-eps. The core of the algorithm 
% is implemented in the function POLEORDERFINDER() below.

%%
% Distance of the sample points from the right end point, i.e. 1.
x = 10.^(-1:-1:-15)';

if ( strcmpi(singEnd, 'right') )       % A pole is expected at x = 1   
    fvalsRight = op(1 - x);
    poleOrder = poleOrderFinder( fvalsRight, x);
else if ( strcmpi(singEnd, 'left') )   % A pole is expected at x = -1
        fvalsLeft = op(-1 + x);
        poleOrder = poleOrderFinder( fvalsLeft, x);
    else
        error('CHEBFUN:SINGFUN:findPoleOrder:unknownPref', ...
              'Blowup preference "%s" unknown.', singEnd )
    end
end

% The algorithm returns a positive number for blow up type singularities. 
% Correct exponents are obtained by negation.
poleOrder = -poleOrder;

end

%%
function poleOrder = poleOrderFinder( fvals, x )
%POLEORDERFINDER   Finds the order of the pole based on function values
%   FVALS given at (1-X).

% We take the absolute value of the function vlaues now and don't call ABS()
% later. This works because scaling with powers of the positive factor (1-x)
% don't change the sign of function values. 
smoothVals = abs(fvals);
if ( any(isinf(smoothVals)) )
    error('CHEBFUN:SINGFUN:findPoleOrder', 'Function returned inf value.')
end

if ( any(isnan(smoothVals)) )
    error('CHEBFUN:SINGFUN:findPoleOrder', 'Function returned NaN value.')
end

% Test parameters
% [TODO]: Should there be a field for this test?
testRatio = 1.01; 
maxPoleOrder = singfun.pref.singfun.maxPoleOrder;

poleOrder = 0;
%smoothVals = smoothVals.*x;

% Loop to see for which power of x the function values become non-divergent
% i.e. when the ratio of function values becomes less then the testRatio.
while ( all(smoothVals(2:end)./smoothVals(1:end-1) > testRatio) && ...
        (poleOrder <= maxPoleOrder) )
    poleOrder = poleOrder + 1;
    smoothVals = smoothVals.*x;
end

if ( poleOrder > maxPoleOrder )
    % Method failed.
    error('CHEBFUN:SINGFUN:fail', ...
            'Pole order excedes limit for maximum pole order.');
end

end
