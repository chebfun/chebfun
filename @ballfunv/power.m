function w = power(v,n)
%.^ Componentwise power for BALLFUNV.
%   F.^G where F is a BALLFUNV and G is a double returns the componentwise 
%   power of F.
%
%   See also BALLFUN/POWER.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check
if isempty( v )
    w = v;
    return
end

V = v.comp;
w = ballfunv( power( V{1},n), power(V{2},n), power(V{3},n) );
end
