function S = size(v)
% SIZE  Size of the expansion coefficients of a BALLFUNV
%   S = SIZE(V) returns a 3 by 3 matrix whose vertical entries represent
%   the sizes of Vx, Vy and Vz
%
%   See also BALLFUN/SIZE.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return the list [m,n,p] corresponding to the size of the ballfun
% functions in f

if isempty( v )
    S = [];
    return
end

V = v.comp;

S = [ size( V{1} );
      size( V{2} );
      size( V{3} ); ];
end
