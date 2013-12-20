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
    %[TODO]: This should be based on tolerances?
    domF = f.funPart.domain;
    domG = g.funPart.domain;
    if ( any( [domF(1) domF(end)] ~= [domG(1) domG(end)] ) ) 
        error( 'CHEBFUN:DELTAFUN:plus', 'f and g must have the same domain' );
    end
    
    s.funPart = f.funPart + g.funPart;
    
    % Add the delta functions:    
    va = f.location;
    vb = g.location;
    A = f.impulses;
    B = g.impulses;
    % Make the number of rows in A and B same by appending zero rows.
    ra = size(A, 1);
    rb = size(B, 1);
    if( ra >= rb )
        B = [ B; zeros(ra-rb, size(B, 2)) ];
    else
        A = [ A; zeros(rb-ra, size(A, 2)) ];
    end      

    s.location = [va, vb];
    s.impulses = [A, B];
    s = simplify(s);       
end

%% %%%%%%%%%%%%% DELTAFUN + FUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive call: If g is a chebfun, upgrade it to a deltafun and call plus
% again.
if ( isa(f, 'deltafun') && isa(g, 'fun') )
    t = deltafun.zeroDeltaFun(g.domain);
    t.funPart = g;
    s = f + t;
    return
end

%% Check for smoothness
% If the result is a smooth object, return a smooth object.
if ( ~anyDelta(s) )
    % return a smooth object if the operation has removed all 
    % delta functions.
    s = s.funPart;
end

end



