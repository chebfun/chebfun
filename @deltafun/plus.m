function s = plus(f, g)
%+   Addition of DELTAFUN objects.
%   F + G adds F and G, where F and G may be DELTAFUN objects or scalars.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
s = deltafun;
% If one of the arguments is empty:
if ( isempty(f) || isempty(g) )
    % return with the empty DELTAFUN:
    return;
end
%%
% One of the arguments i.e. f or g is necessarily a DELTAFUN object. Otherwise, 
% this overloaded plus would not have been called.

% If g is a double, f is a deltafun:
if ( isa(g, 'double') )
    s = f;
    s.funPart = s.funPart + g;
end

% Recursive call, by making the deltafun g as the first argument:
if ( isa(f, 'double') )
    s = g + f;
    return
end

%% DELTAFUN + DELTAFUN
if ( isa(f, 'deltafun') && isa(g, 'deltafun') )
    % If any of f or g have an empty funPart, the funPart of their sum would 
    % be empty, if both are non empty, make sure they have the same domain:
    if( ~isempty(f.funPart) && ~isempty(g.funPart) )
        domF = f.funPart.domain;
        domG = g.funPart.domain;
        if ( any( [domF(1) domF(end)] ~= [domG(1) domG(end)] ) )
            error( 'CHEBFUN:DELTAFUN:plus', 'f and g must have the same domain' );
        end
    end
    s.funPart = f.funPart + g.funPart;
    
    % Add the delta functions:    
    [deltaMag, deltaLoc] = deltafun.mergeImpulses(f.deltaMag, f.location, g.deltaMag, g.location);    
    s.location = deltaLoc;
    s.deltaMag = deltaMag;
end

%% %%%%%%%%%%%%% DELTAFUN + FUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive call: If g is a chebfun, upgrade it to a deltafun and call plus
% again.
if ( isa(f, 'deltafun') && isa(g, 'classicfun') )
    t = deltafun.zeroDeltaFun(g.domain);
    t.funPart = g;
    s = f + t;
    return
end

%% Check for smoothness
% If the result is a smooth object, return a smooth object.
if ( ~anyDelta(s) )
    % Return a smooth object if the operation has removed all 
    % delta functions.
    s = s.funPart;
end

end



