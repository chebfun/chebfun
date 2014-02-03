function f = changeMap(f, newDom)
%MAP   Map the domain of F to NEWDOM via a nonlinear change of variable.
%   F = changeMAP(F, NEWDOM) returns the original F but defined on a new domain
%   [NEWDOM(1) NEWDOM(2)]. The new domain must remain unbounded.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% TODO: Do we require isinf(f.domain) == isinf(newDOM) ?
% TODO: If not, isn't this checked in createMap? 
if ( ~any(isinf(newDom)) )
    error('CHEBFUN:unbndfun:changeMap', ...
        'The new domain must remain unbounded.');
end

% Assign a new nonlinear map to f, obtained from the UNBNDFUN.CREATMAP method:
f.mapping = unbndfun.createMap(newDom);

% Assign the new endpoints of f to the domain property:
f.domain = newDom;

end