function f = plus(f, g)
%+	Addition of two CHEBTECH objects.
%   F + G adds F and G, where F and G may be CHEBTECH objects or scalars.
%
%   If F is a vector-valued chebtech, then F + C is supported if C is a vector
%   of doubles with the same number of columns a F.
%
% See also MINUS, UPLUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % CHEBTECH + [] = []
    
    f = [];

elseif ( isa(g, 'double') ) % CHEBTECH + double
    
    % Update values:
    f.values = bsxfun(@plus, f.values, g); % (For when f and g are vectors).
    % Update coeffs:
    f.coeffs(end,:) = f.coeffs(end,:) + g;
    
    % Update scale:
    vscale = max(f.vscale, max(abs(f.values), [], 1));
    f.epslevel = (f.epslevel*f.vscale + g*eps)./vscale;
    f.epslevel = max(f.epslevel); % [TODO]: Vector epslevel;
    f.vscale = vscale;
    
    
elseif ( isa(f,'double') ) % double + CHEBTECH
    
    % Switch argument order and call CHEBTECH/plus again:
    f = plus(g, f);
    
else % CHEBTECH + CHEBTECH
    
    % Make both CHEBTECH objects have the same length:
    nf1 = size(f.values, 1);
    nf2 = size(g.values, 1);
    if ( nf1 > nf2 )
        % Increase the length of f2 (via prolong):
        g = prolong(g, nf1);
    elseif ( nf1 < nf2 )
        % Increase the length of f1 (via prolong):
        f = prolong(f, nf2);
    end
    
    % Update values and coefficients:
    f.values = f.values + g.values;
    f.coeffs = f.coeffs + g.coeffs;
    
    % Update vscale and epslevel:
    vscale = max(abs(f.values), [], 1);
    f.epslevel = (f.epslevel*f.vscale + g.epslevel*g.vscale)./vscale;
    f.epslevel = max(f.epslevel); % [TODO]: Vector epslevel;
    f.vscale = vscale;
    
    % Look for a zero output:
    if ( ~any(f.values(:)) || ~any(f.coeffs(:)) )
        % Creates a zero CHEBTECH:
        f = f.make(zeros(1, size(f.values, 2)), f.vscale, f.epslevel);
    end
  
end

end
