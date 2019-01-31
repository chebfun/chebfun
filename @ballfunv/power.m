function w = power(v,n)
%.^ Componentwise power for BALLFUNV.
%   F.^G where F is a BALLFUNV and G is a double returns the result from
%   componentwise powers.
%
%   See also BALLFUN/POWER.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

V = v.comp;
w = ballfunv( power( V{1},n), power(V{2},n), power(V{3},n) );
end