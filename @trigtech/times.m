function f = times(f, g, varargin)
%.*   TRIGTECH multiplication.
%   F.*G multiplies TRIGTECH objects F and G or a TRIGTECH by a scalar if either
%   F or G is a scalar.
%
%   If F is an array-valued TRIGTECH, then F.*C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MTIMES, RDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TRIGTECH * [] = []:
if ( isempty(f) || isempty(g) )
    f = []; 
    return
end

if ( ~isa(f, 'trigtech') )      % Ensure F is a TRIGTECH.
    
    f = times(g, f, varargin{:});
    return
    
elseif ( isa(g, 'double') )     % TRIGTECH .* double.
    
    % Do the multiplication:
    if ( size(g, 2) > 1 )
        f.values = bsxfun(@times, f.values, g);
        f.coeffs = f.vals2coeffs(f.values);
    else
        f.values = f.values*g;
        f.coeffs = f.coeffs*g;
    end
    f.isReal = f.isReal & isreal(g);
    return

elseif ( ~isa(f, 'trigtech') || ~isa(g, 'trigtech') )
    % Don't know how to do the operation.
    
    error('CHEBFUN:TRIGTECH:times:typeMismatch', ...
        ['Incompatible operation between objects.\n', ...
         'Make sure functions are of the same type.']);
    
elseif ( size(f.values, 1) == 1 )
    % If we have (constant TRIGTECH).*TRIGTECH, reverse the order and call TIMES
    % again:
    f = times(g, f.values);
    return
    
elseif ( size(g.values, 1) == 1)
    % If we have TRIGTECH.*(constant TRIGTECH), convert the (constant TRIGTECH)
    % to a scalar and call TIMES again:
    f = times(f, g.values); 
    return
end

% Get the size of each TRIGTECH:
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
        error('CHEBFUN:TRIGTECH:times:dim2', ...
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

% Update ishappy:
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
