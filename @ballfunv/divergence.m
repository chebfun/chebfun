function f = divergence(v)
%DIVERGENCE  Divergence of a BALLFUNV in cartesian coordinates.
%   DIVERGENCE(V) is the divergence of the BALLFUNV v expressed in
%   cartesian coordinates.
%
% See also CURL. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( v )
    f = ballfun();
    return
end

[Vx,Vy,Vz] = v.comp{:};
f = diff(Vx,1) + diff(Vy,2) + diff(Vz,3);

end