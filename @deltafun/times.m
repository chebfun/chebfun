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
    h.deltaMag = g * h.deltaMag;
    h = simplify(h);
    return
end

%% Multiplication by a FUN:
if ( isa(g, 'classicfun') )
    % Upgrade to a deltafun and recurse:
    s = deltafun.zeroDeltaFun(g.domain);
    s.funPart = g;
    h = F.*s;
    return
end

%% Multiplication of two DELTAFUNs:
if ( isa(g, 'deltafun') )
    h = deltafun;
    h.funPart = F.funPart .* g.funPart;

    if ( ~isempty(F.location) && ~isempty( g.location) )
       if ( ~isempty(deltafun.numIntersect(F.location, g.location)))
           error( 'CHEBFUN:DELTAFUN:times', 'delta functions at same points can not be multiplied' );
       end
    end
    
    deltaMag1 = [];
    deltaMag2 = [];
    if ( ~isempty(F.location) )
        deltaMag1 = funTimesDelta(g.funPart, F.deltaMag, F.location);
    end
    
    if ( ~isempty(g.location) )
        deltaMag2 = funTimesDelta(F.funPart, g.deltaMag, g.location);
    end
    
    [deltaMag, deltaLoc] = deltafun.mergeImpulses( deltaMag1, F.location, deltaMag2, g.location );
    h.deltaMag = deltaMag;
    h.location = deltaLoc;
else
    % Class of g unknown, throw an error:
    error( 'DELTAFUN:times', 'unknown argument type' );
end
%%
% Check if after multiplication h has become smooth:
if ( ~anyDelta(h) )
    h = h.funPart;
else
    h = simplify(h);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is the distribution f times a  sum of derivatives of delta functions
% all located at the same point x? This is given by the formula:
% < f*delta^(n), phi > = sum < (-1)^(n-j) (n_C_j) * f^(n-j)delta^(j), phi>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = funTimesDelta(f, deltaMag, deltaLoc)

% Make sure there are no redundant rows:
deltaMag = deltafun.cleanRows(deltaMag);
deltaTol = deltafun.pref.deltafun.deltaTol;

% Highest order delta function:
n = size(deltaMag, 1);
m = length(deltaLoc);

% Get all the scaled derivatives needed and store them in a matrix:
Fd = zeros(n, m);
fk = f;
for k = 1:n
    Fd(k, :) = (-1)^(k-1)* feval(fk, deltaLoc);
    fk = diff(fk);
end

D = zeros(n, m);
for j = 1:m
    derivativeCol = Fd(:, j);
    for i = 1:n
        w = zeros(i, 1);
        for k = 1:i
            w(k) = nchoosek(i-1, k-1);
        end
        D(1:i, j) = D(1:i, j) + w.*flipud(derivativeCol(1:i))*deltaMag(i, j);
    end
end

end