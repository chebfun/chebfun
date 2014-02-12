function f = cumsum2( f )
%CUMSUM2   Double indefinite integral of a CHEBFUN2.
%   F = CUMSUM2(F) returns the double indefinite integral of a CHEBFUN2. That is
%                   y  x
%                  /  /
%   CUMSUM2(F) =  |  |   f(x,y) dx dy   for  (x,y) in [a,b] x [c,d],
%                 /  /
%                c  a
%
%   where [a,b] x [c,d] is the domain of f.
% 
% See also CUMSUM, SUM, SUM2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Check for empty:
if ( isempty( f ) ) 
    f = [];
    return
end

f.cols = cumsum( f.cols );   % CUMSUM along the columns.
f.rows = cumsum( f.rows );   % CUMSUM along the rows.

end