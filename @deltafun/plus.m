function s = plus(f, g)
%+   Addition of DELTAFUN objects.
%   F + G adds F and G, where F and G may be DELTAFUN objects or scalars.
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%

% If one of the arguments is empty:
if ( isempty(f) || isempty(g) )
    % Return with an empty DELTAFUN:
    s = deltafun;
    return
end

% One of the arguments (i.e., f or g) is necessarily a DELTAFUN object.
% (Otherwise, this overloaded plus would not have been called!). Ensure it is f:
if ( ~isa(f, 'deltafun') )
    s = plus(g, f);
    return 
end

% Cast g to a DELTAFUN if it is a double or a different kind of FUN:
if ( isa(g, 'double') )
    g = deltafun(0*f.funPart + g);
elseif ( isa(g, 'fun') )
    g = deltafun(g);
end

% Add the funParts:
funPart = f.funPart + g.funPart;

% Add the delta functions:    
[deltaMag, deltaLoc] = deltafun.mergeDeltas(f.deltaMag, f.deltaLoc, ...
    g.deltaMag, g.deltaLoc);    

% Assemble the output DELTAFUN:
s = deltafun(funPart, struct('deltaMag', deltaMag, 'deltaLoc', deltaLoc));

% Simplify:
s = simplifyDeltas(s);

end



