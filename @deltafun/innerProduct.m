function out = innerProduct(f, g)
%INNERPRODUCT Compute the inner product of two DELTAFUN objects.
%   INNERPRODUCT(F, G) is the inner-product of DELTAFUN F and G. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: NH: This requires a test!

%% Trivial cases:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

if ( ~isa(f, 'deltafun') )
    % Ensure first input is a DELTAFUN:
    out = innerProduct(g, f);
    return
end

%% Easy cases:
if ( isa(f, 'deltafun') && isa(g, 'deltafun') ) % Both are DELTAFUNS
    
    % If f and g are both smooth, delegate to the innerproduct of the funPart.
    if ( ~anyDelta(f) && ~anyDelta(g) )   % No deltas
        out = innerProduct(f.funPart, g.funPart);
        return
    elseif ( anyDelta(f) && anyDelta(g) ) % Both deltas
        % If both f and g are non-trivial distributions, then the inner product
        % is not defined. [TODO]: Can we try something in this case as well?
        error('DELTAFUN:innerProduct', ...
            'At least one distribution should be smooth' );
    elseif ( anyDelta(g) )                % Deltas on g
        F = g;
        g = f.funPart;
    else                                  % Deltas on f
        F = f;
    end
    
elseif ( ~anyDelta(f) )                         % f is a (trivial) DELTAFUN
    
    % Both f and g are smooth, so this is the normal innerproduct:
    out = innerProduct(f.funPart, g);
    return
    
else                                            % f is a non-trivial DELTAFUN
    F = f;
end

%% Non-trivial case:
% By this point F is a non-trivial DELTAFUN and g is either a CLASSICFUN or a
% scalar. 
if ( isnumeric(g) )
    if ( isscalar(g) )
        error('DELTAFUN:innerProduct', 'If g is numeric it must be a scalar.');
    end
    % Upgrade a scalar to a CLASSICFUN:
    g = fun.constructor(g, domain(F));
elseif ( ~isa(g, 'fun') )
    error('DELTAFUN:innerProduct', 'Unknown data type.');
end

% Compute the funPart of the inner product.
funIP = innerProduct(F.funPart, g);

% Get location and magnitudes of delta functions:
deltaLoc = F.deltaLoc;
deltaMag = F.deltaMag;
m = size(deltaMag, 1);

% Compute the derivatives needed:
maxDiffOrder = m-1;
G = zeros(m, length(deltaLoc));
G(1,:) = feval(g, deltaLoc);
for k = 1:maxDiffOrder
    g = diff(g);
    G(k+1, :) = feval(g, deltaLoc);
end

% The output is always a scalar double:
deltaIP = 0;
v = ones(maxDiffOrder+1,1); v(2:2:end) = -1;
for k = 1:length(deltaLoc)
    % Apply the definition of inner product with delta functions:
    % <sum(ai dirac^(i)(x-xk)), f(x)> = sum ai*(-1)^i f^(i)(xk)
    ipk = (v.*deltaMag(:, k)).' * G(:, k) ;
    deltaIP = deltaIP + ipk;
end

out = deltaIP + funIP;
end
      
