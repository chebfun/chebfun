function div = divergence(F) 
%DIVERGENCE   Numerical surface divergence of a SPHEREFUNV. 
%   D = DIVERGENCE(F) returns the numerical surface divergence of the
%   SPHEREFUNV. This operations only makes mathematical sense for F that
%   are tanget to the sphere.
%
% See also SPHEREFUNV/DIV, SPHEREFUN/GRAD, SPHEREFUNV/CURL, SPHEREFUNV/VORTICITY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    div = spherefun;
    return
end

% [TODO]: Should we warn the user if F is not tangential to the sphere?
Fc = F.components; 

div = diff(Fc{1}, 1) + diff(Fc{2}, 2) + diff(Fc{3}, 3);

end
