function f = plus(f, g)
%+   Addition of two FOURIERTECH objects.
%   F + G adds F and G, where F and G may be FOURIERTECH objects or scalars.
%
%   If F is an array-valued FOURIERTECH, then F + C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MINUS, UPLUS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % FOURIERTECH + [] = []
    
    f = [];
    
elseif ( isa(g, 'double') ) % FOURIERTECH + double
    
    % Update values (use bsxfun() to handle the case in which g is a vector
    % and f is an array-valued FOURIERTECH):
    f.values = bsxfun(@plus, f.values, g);
    % Update coeffs:
    if ( (size(g, 2) > 1) && (size(f.coeffs, 2) == 1) )
        % Perform singleton expansion of f:
        f.coeffs = repmat(f.coeffs, 1, size(g, 2));
    end
    N = size(f.coeffs,1);
    % Determine the index in the coefficient vector where the constant term
    % is stored. The way we have arranged the coefficients means it should
    % be in the middle of the array, but that depends on whether the number
    % of coefficients is even or odd.
    if mod(N,2) == 1
        const_index = (N+1)/2;
    else
        const_index = N/2;
    end
    f.coeffs(const_index,:) = f.coeffs(const_index,:) + g;
    % Update scale:
    vscaleNew = max(abs(f.values), [], 1);
    % See FOURIERTECH CLASSDEF file for documentation on this:
    f.epslevel = (f.epslevel.*f.vscale + abs(g)*eps)./vscaleNew;
    f.vscale = vscaleNew;
    
elseif ( isa(f, 'double') ) % double + FOURIERTECH
    
    % Switch argument order and call FOURIERTECH/PLUS again:
    f = plus(g, f);
    
else % FOURIERTECH + FOURIERTECH
    
    % We will simply add the values together then compute the coefficients
    % of the result.  This is probably not the most efficient means of
    % determing the sum.
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
    f.coeffs = f.vals2coeffs(f.values);
    
    % Look for a zero output:
    if ( ~any(f.values(:)) || ~any(f.coeffs(:)) )
        % Create a zero FOURIERTECH:
        epslevel = max(f.epslevel, g.epslevel);
        ishappy = f.ishappy && g.ishappy;
        z = zeros(1, size(f.values, 2));
        f = f.make(z, z, f.hscale);
        f.epslevel = epslevel;
        f.ishappy = ishappy;
    else
        % Update vscale, epslevel, and ishappy:
        vscaleNew = max(abs(f.values), [], 1);
        % See FOURIERTECH CLASSDEF file for documentation on this:
        f.epslevel = (f.epslevel.*f.vscale + g.epslevel.*g.vscale)./vscaleNew;
        f.vscale = vscaleNew;
        f.ishappy = f.ishappy && g.ishappy;
    end
    
end

end
