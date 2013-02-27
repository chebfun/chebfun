function f = times(f, g, varargin)
%.*	FUNCHEB multiplication.
%   F.*G multiplies FUNCHEB objects F and G or a FUNCHEB by a scalar
%   if either F or G is a scalar.
%
% See also MTIMES, RDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% FUNCHEB * [] = []
if ( isempty(f) || isempty(g) )
    f = []; 
    return
end

if ( ~isa(f, 'funcheb') )      % Ensure F is a FUNCHEB
    
    f = times(g, f, varargin{:});
    return
    
elseif ( isa(g, 'double') )     % FUNCHEB * double
    
    % Check dimensions:
    if ( (size(g, 1) > 1) || ...
         ((size(g, 2) > 1) && (size(f.values, 2) ~= size(g, 2))) )
        error('CHEBFUN:FUNCHEB:times:dim', ...
            'Inner matrix dimensions must agree.');
    end
    
    % Do the multiply:
    if ( size(g, 2) > 1 )
        f.values = bsxfun(@times, f.values, g);
        f.coeffs = bsxfun(@times, f.coeffs, g);
        f.vscale = f.vscale.*abs(g);
    else
        f.values = f.values*g;
        f.coeffs = f.coeffs*g;
        f.vscale = f.vscale*abs(g);
    end
    
    return
    
elseif ( size(f.values, 1) == 1 ) % constant FUNCHEB * FUNCHEB
    f = times(g, f.values);
    f.epslevel = max(f.epslevel, g.epslevel);
    return
    
elseif ( size(g.values, 1) == 1)  % FUNCHEB * constant FUNCHEB
    f = times(f, g.values); 
    f.epslevel = max(f.epslevel, g.epslevel);
    return
    
end

% Determine a tolerance if none is given:
if ( nargin < 3 )
    pref = f.pref; 
else
    pref = varargin{1};
end

% Get the size of each funcheb:
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
        error('CHEBFUN:FUNCHEB:times:dim2', ...
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
% vscale = norm([f.vscale ; g.vscale ; norm(values(:), inf)], inf);
vscale = max(f.vscale, g.vscale);
vscale = max(vscale, max(abs(values)));
epslevel = max(f.epslevel, g.epslevel);

% Assign back to f:
f.values = values;
f.coeffs = f.chebpoly(values);
f.vscale  = vscale;
f.epslevel = epslevel;

% Simplify!
f = simplify(f, pref);

% Funs F and G are such that their product should be positive. Enforce this on
% the values. (Simplify could have ruined this property).
if ( pos )
    f.values = abs(f.values); 
    f.coeffs = f.chebpoly(f.values);
end

end
