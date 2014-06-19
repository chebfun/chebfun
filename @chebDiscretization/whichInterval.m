function intNum = whichInterval(disc, location, direction)
%WHICHINTERVAL Index of the subinterval of a given point and direction. 
%   INTNUM = WHICHINTERVAL(DISC, LOC, DIRN) returns which interval within
%   DISC.domain a given point LOC is located in. The DIRN is +1 or -1 and
%   indicates approach from the right or left (zero means don't care).

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: How does this differ from chebfun.whichInterval?

% Quick shortcut if there are no breakpoints.
numIntervals = length(disc.domain);
if (  numIntervals == 2 )
    intNum = 1;
    return
end

% If there were no direction, this would be it. 
intNum = find( location >= disc.domain, 1, 'last' );

% LINOP already screened to make sure the location is in the interval, so we can
% check for being at the right endpoint.
if ( intNum == numIntervals )
    if ( direction > 0 )
        error('CHEBFUN:CHEBDISCRETIZATION:whichInterval:undefined', ...
            'Evaluation direction is undefined at the location.')
    end
    direction = -1;  % this forces the adjustment below
end

% Need to decrement if at a breakpoint coming from the left, or if at
% the right endpoint.
len = disc.domain(end) - disc.domain(1);
if ( direction < 0 && abs( location - disc.domain(intNum) ) < 10*eps*len )
    intNum = intNum - 1;
end

end
