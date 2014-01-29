function out = innerProduct(f, g)
%INNERPRODCUT Compute the inner product of two DELTAFUN objects.
%   INNERPROTCUT(F, G) is the inner-product of DELTAFUN F and G. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Empty cases:
if( isempty(f) || isempty(g) )
    out = [];
    return
end

if( isa(f, 'deltafun') && isa(g, 'deltafun') )
    % If f and g are both smooth, delegate to
    % the innerproduct of the function part.
    if( ~anyDelta(f) && ~anyDelta(g) )
        out = f.funPart' * g.funPart;
        return
    end
    % If both f and g are non-trivial distributions, then the inner product is
    % not defined.
    % [TODO]: Can we try something in this case as well?
    if ( anyDelta(f) && anyDelta(g) )
        error('DELTAFUN:innerProduct', 'At least one distribution should be smooth' );
    end    
end

% One deltafun and the other something else. Make sure F is a deltafun and 
% g contains the other argument:
if( isa(f, 'deltafun') )
    F = f;    
else
    % If f is not a DELTAFUN, g must be:
    F = g;
    g = f;
end

% If F is smooth, then since g is also smooth, this is the nomral innerproduct.
if ( ~anyDelta(F) )
    out = F.funPart' * g;
    return;
end

% By this point, F is not a smooth DELTAFUN. If g is not a CLASSICFUN or 
% a DELTAFUN, we upgrade it to a DELTAFUN with constant FUNPART, if g is
% something else, it's an error.
if( isnumeric(g) )
    if( numel(g) > 1 )
        error( 'DELTAFUN:innerProduct', 'if g is numeric, it should be a scalar' );
    end
    g = fun.constructor(g, F.funPart.domain);
elseif( ~isa(g, 'fun') )
    error( 'DELTAFUN:innerProduct', 'unknown data type' );
end

% Compute the smooth part of the inner product.
if ( ~isempty(F.funPart) )
    funIP = innerProduct(F.funPart, g);
else
    funIP = 0;
end

% Get location and magnitudes of delta functions:
f = simplify(f);
deltaLoc = f.deltaLoc;
numDeltas = length(deltaLoc);
deltaMag = f.deltaMag;
m = size(deltaMag, 1);

% Compute the derivatives needed:
maxDiffOrder = m-1;

G = zeros(m, numDeltas);
G(1,:) = feval(g, deltaLoc);
for k = 1:maxDiffOrder
    g = diff(g);
    G(k+1, :) = feval(g, deltaLoc);
end

% The output is always a scalar double:
deltaIP = 0;
v = (-1).^(0:maxDiffOrder).';
for k = 1:length(deltaLoc)
    %[TODO]: Please explain this voodoo
    ipk = (v.*deltaMag(:, k)).' * G(:, k) ;
    deltaIP = deltaIP + ipk;
end

out = deltaIP + funIP;
end
      