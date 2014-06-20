function f = times(f, g, varargin)
%.*   FOURTECH multiplication.
%   F.*G multiplies FOURTECH objects F and G or a FOURTECH by a scalar if either
%   F or G is a scalar.
%
%   If F is an array-valued FOURTECH, then F.*C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MTIMES, RDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% FOURTECH * [] = []:
if ( isempty(f) || isempty(g) )
    f = []; 
    return
end

if ( ~isa(f, 'fourtech') )      % Ensure F is a FOURTECH.
    
    f = times(g, f, varargin{:});
    return
    
elseif ( isa(g, 'double') )     % FOURTECH .* double.
    
    % Do the multiplication:
    if ( size(g, 2) > 1 )
        f.values = bsxfun(@times, f.values, g);
        f.coeffs = f.vals2coeffs(f.values);
        f.vscale = f.vscale.*abs(g);
    else
        f.values = f.values*g;
        f.coeffs = f.coeffs*g;
        f.vscale = f.vscale*abs(g);
    end
    f.epslevel = f.epslevel + eps(g);
    f.isReal = f.isReal & isreal(g);
    return

elseif ( ~isa(f, 'fourtech') || ~isa(g, 'fourtech') )
    % Don't know how to do the operation.
    
    error('CHEBFUN:FOURTECH:times:typeMismatch', ...
        ['Incompatible operation between objects.\n', ...
         'Make sure functions are of the same type.']);
    
elseif ( size(f.values, 1) == 1 )
    % If we have (constant FOURTECH).*FOURTECH, reverse the order and call TIMES
    % again:
    f = times(g, f.values);
    f.epslevel = max(f.epslevel, g.epslevel);
    return
    
elseif ( size(g.values, 1) == 1)
    % If we have FOURTECH.*(constant FOURTECH), convert the (constant FOURTECH)
    % to a scalar and call TIMES again:
    f = times(f, g.values); 
    f.epslevel = max(f.epslevel, g.epslevel);
    return
end

% Get the size of each FOURTECH:
[fn, fm] = size(f.values);
[gn, gm] = size(g.values);

% The length of the product is known:
fNew = prolong(f, fn + gn - 1); 

% Check dimensions:
if ( fm ~= gm )
    if ( fm == 1 )
        % Allow [Inf x 1] .* [Inf x m].
        fNew.values = repmat(fNew.values, 1, gm);
        fNew.coeffs = repmat(fNew.coeffs, 1, gm);
    elseif ( gm == 1 )
        % Allow [Inf x m] .* [Inf x 1].
        g.values = repmat(g.values, 1, fm);
        g.coeffs = repmat(g.coeffs, 1, fm);
    else
        error('CHEBFUN:FOURTECH:times:dim2', ...
            'Inner matrix dimensions must agree.');
    end
end

% Check for two cases where the output is known in advance to be positive,
% namely F == conj(G) or F == G and isreal(F).
pos = false;

% Multiply values:
if ( isequal(f, g) )
   values = fNew.values.^2;          
   if ( isreal(f) )
       pos = true; 
   end
elseif ( isequal(conj(f), g) )
   values = conj(fNew.values).*fNew.values;
   pos = true;
else
   gNew = prolong(g, fn + gn - 1); 
   values = fNew.values.*gNew.values;
end

% Assign values and coefficients back to f:
f.values = values;
f.coeffs = f.vals2coeffs(values);

% Update vscale, epslevel, and ishappy:
vscale = max(abs(f.values), [], 1);

% Avoid NaNs:
tmpVscale = vscale;
tmpVscale(vscale == 0) = 1;
f.vscale(f.vscale == 0) = 1;
g.vscale(g.vscale == 0) = 1;
 
% See CHEBTECH CLASSDEF file for documentation on this:
epslevelBound = (f.epslevel + g.epslevel) .* (f.vscale.*g.vscale./tmpVscale);
f.epslevel = updateEpslevel(f, epslevelBound);
f.vscale  = vscale;
f.ishappy = f.ishappy && g.ishappy;

% Simplify.
f = simplify(f);

if ( pos )
    % Here we know that the product of F and G should be positive. However,
    % SIMPLIFY may have destroyed this property, so we enforce it.
    f.values = abs(f.values); 
    f.coeffs = f.vals2coeffs(f.values);
    f.isReal = true(1, size(f.coeffs, 2));
else
    f.isReal = f.isReal & g.isReal;
end

% If f and g are real then make the result real.
f.values(:,f.isReal) = real(f.values(:,f.isReal));

end