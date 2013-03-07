function f = plus(f, g)
%+	Addition of two FUNCHEB objects.
%   F + G adds F and G, where F and G may be FUNCHEB objects or scalars.
%
%   If F is a vector-valued funcheb, then F + C is supported if C is a vector
%   of doubles with the same number of columns a F.
%
% See also MINUS, UPLUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % FUNCHEB + [] = []
    
    f = [];

elseif ( isa(g, 'double') ) % FUNCHEB + double
    
    % Update values:
    f.values = bsxfun(@plus, f.values, g); % (For when f and g are vectors).
    % Update coeffs:
    f.coeffs(end,:) = f.coeffs(end,:) + g;
    % Update scale:
    f.vscale = max(f.vscale, max(abs(f.values), [], 1));
    
elseif ( isa(f,'double') ) % double + FUNCHEB
    
    % Switch argument order and call FUNCHEB/plus again:
    f = plus(g, f);
    
else % FUNCHEB + FUNCHEB
    
    % Make both FUNCHEB objects have the same length:
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
    
    % Update epslevel:
    f.epslevel = max(f.epslevel, g.epslevel);
    
    % Update scales:
    f.vscale = max(max(f.vscale, g.vscale), max(abs(f.values), [], 1));
    
    % Look for a zero output:
    if ( ~any(f.values(:)) || ~any(f.coeffs(:)) )
        % Creates a zero FUNCHEB:
        f = f.make(zeros(1, size(f.values, 2)), f.vscale, f.epslevel);
    end
  
end

end
