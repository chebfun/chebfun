function f = plus(f, g)
%+   Addition of two CHEBTECH objects.
%   F + G adds F and G, where F and G may be CHEBTECH objects or scalars.
%
%   If F is an array-valued CHEBTECH, then F + C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MINUS, UPLUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % CHEBTECH + [] = []
    
    f = [];
    
elseif ( isa(g, 'double') ) % CHEBTECH + double
    
    % Update values (use bsxfun() to handle the case in which g is a vector
    % and f is an array-valued CHEBTECH):
    f.values = bsxfun(@plus, f.values, g);
    % Update coeffs:
    if ( (size(g, 2) > 1) && (size(f.coeffs, 2) == 1) )
        % Perform singleton expansion of f:
        f.coeffs = repmat(f.coeffs, 1, size(g, 2));
    end
    f.coeffs(end,:) = f.coeffs(end,:) + g;
    % Update scale:
    vscale = max(abs(f.values), [], 1);
    % See CHEBTECH CLASSDEF file for documentation on this:
    f.epslevel = (f.epslevel*f.vscale + abs(g)*eps)./vscale;
    f.epslevel = max(f.epslevel); % [TODO]: Vector epslevel;
    f.vscale = vscale;
    
elseif ( isa(f, 'double') ) % double + CHEBTECH
    
    % Switch argument order and call CHEBTECH/PLUS again:
    f = plus(g, f);
    
else % CHEBTECH + CHEBTECH
    
    % Make both CHEBTECH objects have the same length:
    nf = size(f.values, 1);
    ng = size(g.values, 1);
    if ( nf > ng )
        % Increase the length of g (via PROLONG):
        g = prolong(g, nf);
    elseif ( nf < ng )
        % Increase the length of f (via PROLONG):
        f = prolong(f, ng);
    end
    
    % Update values and coefficients:
    f.values = f.values + g.values;
    f.coeffs = f.coeffs + g.coeffs;
    
    % Look for a zero output:
    if ( ~any(f.values(:)) || ~any(f.coeffs(:)) )
        % Create a zero CHEBTECH:
        epslevel = max(f.epslevel, g.epslevel);
        ishappy = f.ishappy && g.ishappy;
        z = zeros(1, size(f.values, 2));
        f = f.make(z, z, f.hscale);
        f.epslevel = epslevel;
        f.ishappy = ishappy;
    else
        % Update vscale, epslevel, and ishappy:
        vscale = max(abs(f.values), [], 1);
        % See CHEBTECH CLASSDEF file for documentation on this:
        f.epslevel = (f.epslevel*f.vscale + g.epslevel*g.vscale)./vscale;
        f.epslevel = max(f.epslevel);  % [TODO]: Vector epslevel;
        f.vscale = vscale;
        f.ishappy = f.ishappy && g.ishappy;
    end
    
end

end
