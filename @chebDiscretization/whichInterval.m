function intnum = whichInterval(disc,location,direction)
%WHICHINTERVAL Index of the subinterval of a given point and direction. 
%   WHICHINTERVAL(DISC,LOC,DIRN) finds which interval within DISC.domain a given
%   point LOC is located in. The DIRN is +1 or -1 and indicates approach from
%   the right or left (zero means don't care). 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Quick shortcut if there are no breakpoints.
numIntervals = length(disc.domain);
if (  numIntervals == 2 )
    intnum = 1;
    return
end

% If there were no direction, this would be it. 
intnum = find( location >= disc.domain, 1, 'last' );

% linop already screened to make sure the location is in the
% interval, so we can check for being at the right endpoint.
if ( intnum == numIntervals )
    if (direction > 0)
        error('Evaluation direction is undefined at the location.')
    end
    direction = -1;  % this forces the adjustment below
end

% Need to decrement if at a breakpoint coming from the left, or if at
% the right endpoint.
len = disc.domain(end) - disc.domain(1);
if ( direction < 0 && abs( location - disc.domain(intnum) ) < 10*eps*len )
    intnum = intnum - 1;
end

end
