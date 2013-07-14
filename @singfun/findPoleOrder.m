function poleOrder = findPoleOrder(op, singEnd)
% Given a function f defined on [-1,1], find an integer exponent
% E such that f*(1-x)^(-E) is bounded at 1-eps.

% distance of sample points from the end points
x = 10.^(-1:-1:-15)';
% if a pole is expected at x = 1
if ( strcmpi(singEnd, 'right') )
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

end

function poleOrder = poleOrderFinder( fvals, x )
% Iteratively increase the proposed value. If
% the singularity at the right endpoint can be mollified, we can declare
% success.

% we take the absolute value here, 
% and since (1-x) is assumed positive,
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
    warning('CHEBFUN:SINGFUN:fail',...
        'Could not detect function singularity.');
    poleOrder = 0;
end
end

