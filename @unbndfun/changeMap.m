function f = changeMap(f, newdom)
%MAP   Maps the domain of an UNBNDFUN F to NEWDOM via a nonlinear change of 
%      variable.
%   F = changeMAP(F,NEWDOM), where the UNBNDFUN F has a semi-infinite or 
%   infinite domain, returns the original F but defined on a new domain 
%   [NEWDOM(1) NEWDOM(2)].

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( ~any(isinf(newdom)) )
    error('CHEBFUN:unbndfun:changeMap', 'The new domain must remain unbounded.');
end

% Assign a new nonlinear map to f, obtained from the UNBNDFUN.CREATMAP method:
f.mapping = unbndfun.createMap(newdom);

% Assign the new endpoints of f to the domain property:
f.domain = newdom;

end