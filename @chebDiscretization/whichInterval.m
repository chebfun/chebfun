function intnum = whichInterval(disc,location,direction)

% Which subinterval of the domain is active?
intnum = find( location >= disc.domain, 1, 'last' );

% linop already screened to make sure the location is in the
% interval, so we can check for being at the right endpoint.
if ( intnum == length(disc.domain) )
    if (direction > 0)
        error('Evaluation direction is undefined at the location.')
    end
    direction = -1;  % this forces the adjustment below
end

% Need to adjust if at a breakpoint coming from the left, or if at
% the right endpoint.
len = disc.domain(end) - disc.domain(1);
if ( direction < 0 && abs( location - disc.domain(intnum) ) < 10*eps*len )
    intnum = intnum - 1;
end

end
