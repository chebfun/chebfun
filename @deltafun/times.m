function h = times(f, g)
%.*   Multiply DELTAFUNS with DELTAFUNS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: This method will be called only if both F and G are DELTAFUNS or at the
% most one of F and G is a scalar double.

%% Empty case:
h = deltafun;
if ( isempty(f) || isempty(g) )
    return
end

% Make sure F is a DELTAFUN and copy the other input in g:
if ( ~isa(f, 'deltafun') )
    h = times(g, f);
    return
end

%% Multiplication by a scalar
if ( isa(g, 'double') )
    h = f;
    % Scale everything and return:
    h.funPart = g * f.funPart;
    h.deltaMag = g * h.deltaMag;
    h = simplifyDeltas(h);
    return
end

%% Multiplication by a FUN:
if ( isa(g, 'classicfun') )
    % Upgrade to a DELTAFUN:
    g = deltafun(g, [], []);
end

%% Multiplication of two DELTAFUNs:
if ( isa(g, 'deltafun') )
    funPart = f.funPart .* g.funPart;
    
    if ( ~isempty(f.deltaLoc) && ~isempty( g.deltaLoc) )
        if ( ~isempty(deltafun.numIntersect(f.deltaLoc, g.deltaLoc)))
            error( 'CHEBFUN:DELTAFUN:times:times', ...
                'Delta functions at same points can not be multiplied' );
        end
    end
    
    if ( ~isempty(f.deltaLoc) )
        deltaMag1 = funTimesDelta(g.funPart, f.deltaMag, f.deltaLoc);
    else
        deltaMag1 = [];
    end
    if ( ~isempty(g.deltaLoc) )
        deltaMag2 = funTimesDelta(f.funPart, g.deltaMag, g.deltaLoc);
    else
        deltaMag2 = [];
    end
    
    % Merge delta functions:
    [deltaMag, deltaLoc] = ...
        deltafun.mergeDeltas(deltaMag1, f.deltaLoc, deltaMag2, g.deltaLoc);
    
    % Assumble the DELTAFUN:
    h = deltafun(funPart, struct('deltaMag', deltaMag, 'deltaLoc', deltaLoc));
else
    % Class of g unknown, throw an error:
    error( 'CHEBFUN:DELTAFUN:times:unknownType', 'unknown argument type' );
end

%%
% Check if after multiplication h has become smooth:
h = simplifyDeltas(h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is the distribution f times a  sum of derivatives of delta functions
% all located at the same point x? This is given by the formula:
% < f*delta^(n), phi > = sum < (-1)^(n-j) (n_C_j) * f^(n-j)delta^(j), phi>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = funTimesDelta(f, deltaMag, deltaLoc)

% Make sure there are no redundant rows:
deltaMag = deltafun.cleanRows(deltaMag);

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
