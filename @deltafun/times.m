function h = times(f,g)
%.*   Multiply DELTAFUNS with DELTAFUNS

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: This method will be called only if both F and G are DELTAFUNS or at the
% most one of F and G is a scalar double.

%% Empty case:
h = deltafun;
if ( isempty(f) || isempty(g) )
    return
end

% Trivial cases:
% Check if inputs are other than DELTAFUNS or doubles:
if ( (~isa(f, 'deltafun') && ~isa(f, 'double')) || ...
     (~isa(g, 'deltafun') && ~isa(g, 'double')) )
    error( 'DELTAFUN:times', 'Input can only be a DELTAFUN or a double' )
end

% Make sure F is a deltafun and copy the other input in g
if ( ~isa(f, 'deltafun') )
    % Then g must be a deltafun
    F = g;
    g = f;
else
    % f is a deltafun
    F = f;
end
    
%% Multiplication by a scalar
if ( isa(g, 'double') )
    h = F;
    % Scale everything and return:
    h.funPart = g * F.funPart;
    h.impulses = g * h.impulses;
    return
end

%% Multiplication by a CHEBFUN:
if ( isa(g, 'chebfun') )
    % Upgrade to a deltafun and recurse:
    s = deltafun.zeroDeltaFun(g.domain);
    s.funPart = g;
    h = F.*s;
    return
end

%% Multiplication of two DELTAFUNs:
if ( isa(g, 'deltafun') )
    h = F;
    h.funPart = F.funPart .* g.funPart;
    if ( ~isempty( g.location) )
       if ( ~isempty(deltafun.numIntersect(F.location, g.location)))
           error( 'CHEBFUN:DELTAFUN:times', 'delta functions at same points can not be multiplied' );
       else
           % The magic goes here:
           
       end
       
    end
else
    % Class of g unknown, throw an error
    error( 'CHEBFUN:DELTAFUN:times', 'unknown argument type' );
end
%%
% Check if after multiplication h has become smooth:
if ( issmooth(h) )
    h = h.smoothPart;
end

end
