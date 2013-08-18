function poleOrder = findPoleOrder(op, singEnd)
%FINDPOLEORDER   Finds the order of the pole in the function handle OP at
%   x = 1 or -1 depending upon the string 'left' or 'right' passed in 
%   SINGEND.
%   
% Example:
%   p = singfun.findPoleOrder(@(x) 1./(1-x), 'right' )
%   p = singfun.findPoleOrder(@(x) 1./(1+x).^2, 'left' )
%
% See also SINGFUN.FINDSINGORDER and SINGFUN.FINDSINGEXPONENTS.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% The algorithm tries to find the smallest non-negative integer
% k such that op(x).*(1-x)^k is bounded at 1-eps. The algorithm is actually
% implemented in the function POLEORDERFINDER() below.

%%
% distance of the sample points from the right end point, i.e. 1.
x = 10.^(-1:-1:-15)';

if ( strcmpi(singEnd, 'right') )
    % if a pole is expected at x = 1
    fvalsRight = op(1-x);
    poleOrder = poleOrderFinder( fvalsRight, x);
else if ( strcmpi(singEnd, 'left') )
        % if a pole is expected at x = -1
        fvalsLeft = op(-1+x);
        poleOrder = poleOrderFinder( fvalsLeft, x);
    else
        error('CHEBFUN:SINGFUN:findPoleOrder:unknownPref',...
                    'Blowup preference "%s" unknown', singEnd )
    end
end

% return negative numbers as exponents
poleOrder = -poleOrder;
end
%%
function poleOrder = poleOrderFinder( fvals, x )
%POLEORDERFINDER   Finds the order of the pole based on function values
%   FVALS given at (1-X).

% we take the absolute value of the function vlaues now and since the 
% factor (1-x) is assumed positive, we dont' need ABS() afterwards when 
% function values are scaled by powers of (1-x).
smoothVals = abs(fvals);
if ( any(isinf(smoothVals)) )
    error( 'CHEBFUN:SINGFUN:findPoleOrder', 'function returned inf while evaluation')
end

if ( any(isnan(smoothVals)) )
    error( 'CHEBFUN:SINGFUN:findPoleOrder', 'function returned NaN while evaluation')
end

% Test parameters
testRatio = 1.01;
maxPoleOrder = singfun.pref.singfun.maxPoleOrder;

poleOrder = 0;
%smoothVals = smoothVals.*x;
% Loop to see for which power of x the function values become non-divergent
% i.e. when the ratio of function values becomes less then the testRatio.
while( all(smoothVals(2:end)./smoothVals(1:end-1) > testRatio) && ( poleOrder <= maxPoleOrder ) )
    poleOrder = poleOrder + 1;
    smoothVals = smoothVals.*x;
end
if ( poleOrder > maxPoleOrder )
    % Method failed.
    % [TODO]: Error may be?
    warning('CHEBFUN:SINGFUN:fail',...
        'Pole order excedes limit for maximum pole order.');
    poleOrder = 0;
end

end