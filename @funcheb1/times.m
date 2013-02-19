function f = times(f, g, varargin)
%.*	FUNCHEB1 multiplication.
%   F.*G multiplies FUNCHEB1 objects F and G or a FUNCHEB1 by a scalar
%   if either F or G is a scalar.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% FUNCHEB1 * [] = []
if ( isempty(f) || isempty(g) )
    f = []; 
    return
end

if ~isa(f, 'funcheb1')          % Ensure F is a FUNCHEB1
    
    f = times(g, f, varargin{:});
    return
    
elseif isa(g, 'double')         % FUNCHEB1 * double
    
    % Check dimensions:
    if ( size(g, 1) > 1 || size(f.values, 2) ~= size(g, 2) )
        error('CHEBFUN:FUNCHEB1:times:dim', ...
            'Inner matrix dimensions must agree.');
    end
    
    % Do the multiply:
    if size(g, 2) > 1
        n = length(f.values);
        f.values = f.values.*repmat(g, n, 1);
        f.coeffs = f.coeffs.*repmat(g, n, 1);
        f.vscale = f.vscale*norm(g, inf);
    else
        f.values = f.values*g;
        f.coeffs = f.coeffs*g;
        f.vscale = f.vscale*abs(g);
    end
    
    return
    
elseif size(f.values, 1) == 1   % constant FUNCHEB1 * FUNCHEB1
    
    f = mtimes(f.values, g);
    f.epslevel = max(f.epslevel, g.epslevel);
    return
    
elseif size(g.values,1 ) == 1   % FUNCHEB1 * constant FUNCHEB1
    
    f = mtimes(f, g.values); 
    f.epslevel = max(f.epslevel, g.epslevel);
    return
    
end

% Determine a tolerance if none is given:
if ( nargin < 3 )
    pref = funcheb1.pref; 
else
    pref = varargin{1};
end

% Get the length of each funcheb1:
[fn, fm] = size(f.values);
[gn, gm] = size(g.values);

% The length of the product is known:
fNew = prolong(f, fn + gn - 1); 

% Check dimensions:
if ( fm ~= gm )
    if ( fm == 1 )
        % Allow [inf x 1] .* [inf x m].
        fNew.values = repmat(fNew.values, 1, gm);
        fNew.coeffs = repmat(fNew.coeffs, 1, gm);
    elseif ( gm == 1 )
        % Allow [inf x m] .* [inf x 1].
        g.values = repmat(g.values, 1, fm);
        g.coeffs = repmat(g.coeffs, 1, fm);
    else
        error('CHEBFUN:FUNCHEB1:times:dim2', ...
            'Inner matrix dimensions must agree.');
    end
end

% Determine if output should be positive
% (i.e., F == conj(G) or F == G and isreal(F)).
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

% Update scales:
vscale = norm([f.vscale ; g.vscale ; norm(values(:), inf)], inf);
epslevel = max(f.epslevel, g.epslevel);

% Simplify (and obtain coeffs):
pref.funcheb1.eps = max(pref.funcheb1.eps, epslevel);
coeffs = funcheb1.chebpoly(values);
[cutoff, epslevel] = funcheb1.checkCoeffs(values, coeffs, vscale, pref); 
coeffs = funcheb1.alias(coeffs, cutoff);
values = funcheb1.chebpolyval(coeffs);

% Funs F and G are such that their product should be positive.
% Enforce this on the values. (Simplify could have ruined this property).
if ( pos )
    values = abs(values); 
end

% Assign values to FUNCHEB2:
f.values = values;
f.coeffs = coeffs;
f.epslevel = epslevel;
