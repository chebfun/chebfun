function poleOrder = findPoleOrder(op, singFlag)
% finds the order of the poles of the operator OP at x = -1 and x = 1.
% SINGFLAG is a 2x1 vector. A NAN indicates whehter a pole at the the point x =-1
% or x =1 is expected.
% Given a function f defined on [-1,1], find an integer exponent
% E such that f*(1-x)^(-E) is bounded at 1-eps.

poleOrder = zeros(1,2);

% distance of sample points from the end points
x = 10.^(-1:-1:-15)';
% if a pole is expected at x = 1
if ( isnan( singFlag(2) ) )
    fvalsRight = op(1-x);
    poleOrder(2) = poleOrderFinder( fvalsRight, x );
else
    % no singularity
    poleOrder(2) = 0;
end
% if a pole is expected at x = -1
if ( isnan( singFlag(1) ) )  
    fvalsLeft = op(-1+x);
    poleOrder(1) = poleOrderFinder( fvalsLeft, x );
else
    % no singularity
    poleOrder(1) = 0;
end
end

function poleOrder = poleOrderFinder( fvals, x )
% Iteratively increase the proposed value. If
% the singularity at the right endpoint can be mollified, we can declare
% success.

% we take the absolute value here, 
% and since x is assumed positive,
% we dont' need abs() afterwards.
smoothVals = abs(fvals);
poleOrder = 0;

% Test parameters
testRatio = 1.01;
maxPoleOrder = 100;

% Loop to see for which power of x the
% function values become non-divergent
while( all(smoothVals(2:end)./smoothVals(1:end-1) > testRatio ) && ( poleOrder <= maxPoleOrder ) )
    poleOrder = poleOrder + 1;
    smoothVals = smoothVals.*x;
end
if ( poleOrder > maxPoleOrder )
    % Failure mode: do nothing.
    warning('CHEBFUN:singfun:fail',...
        'Could not detect function singularity.');
    poleOrder = 0;
end
end