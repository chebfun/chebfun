function out = ip(f, g)
%IP   Inner-product of two distributions.
%   IP(F, G) At least one of F or G has to be smooth or a chebfun.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Empty cases:
if( isempty(f) || isempty(g) )
    out = [];
    return
end

if( isa(f, 'deltafun') && isa(g, 'deltafun') )
    % If f and g are both smooth, delegate to
    % innerproduct of the function part.
    if( issmooth(f) && issmooth(g) )
        out = f.funPart' * g.funPart;
        return
    end
    % If both f and g are non-trivial distributions, then the inner product is
    % not defined.
    % [TODO]: We can try something in this case as well.
    if ( ~issmooth(f) && ~issmooth(g) )
        error('CHEBFUN:DELTAFUN:ip', 'At least one distribution should be smooth' );
    end    
end

% One deltafun and the other something else. Make sure F is a deltafun and 
% g contains the other argument:
if( isa(f, 'deltafun') )
    F = f;    
else
    % g must be a deltafun
    F = g;
    g = f;
end

% If F is smooth, then since g is also smooth, this is the nomral innerproduct.
if ( issmooth(F) )
    out = F.funPart' * g;
    return;
end

% By this point, F is not smooth and g is something. If g is a chebfun, we don't
% do anything, if g is a double, we upgrade it to a constant chebfun, if g is
% something else, it is an error.
if( isnumeric(g) )
    if( numel(g) > 1 )
        error( 'CHEBFUN:DELTAFUN:ip', 'if g is numeric, it should be a scalar' );
    end
    g = chebfun(g, F.funPart.domain);
elseif( ~isa(g, 'chebfun') )
    error( 'CHEBFUN:DELTAFUN:ip', 'unknown data type' );
end


% Compute the smooth part of the inner product.
smoothIP = F.funPart'*g;

% Computing the inner product

% Compute the derivatives needed:
diffOrder = F.delta.diffOrder;
maxDiffOrder = max(diffOrder);
G = cell(maxDiffOrder, 1);
G{1} = g;
for k = 1:maxDiffOrder
    g = diff(g);
    G{k+1} = g;
end

% The output is always a scalar double
deltaIP = 0;

% Get location and magnitudes of delta functions:
deltaLoc = f.delta.location;
deltaMag = f.delta.magnitude;

for k = 1:length(deltaLoc)    
    gk = (-1)^diffOrder(k) * G{diffOrder(k)+1};
    gkVal = deltaMag(k) * gk(deltaLoc(k));
    
    if( f.delta.isReal(k) )
        gkVal = real(gkVal);
    end
    
    if( f.delta.isImag(k) )
        gkVal = imag(gkVal);
    end
    
    if( f.delta.isConj(k) )
        gkVal = conj(gkVal);
    end
    deltaIP = deltaIP + gkVal;
end

out = deltaIP + smoothIP;
end