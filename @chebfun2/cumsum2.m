function f = cumsum2(f)
%CUMSUM2 Double indefinite integral of a chebfun2.
%
% F = CUMSUM2(F) returns the double indefinite integral of a chebfun2. That
% is,
%                  y  x
%                 /  /
%  CUMSUM2(F) =  |  |   f(x,y) dx dy   for  (x,y) in [a,b] x [c,d],
%                /  /
%               c  a
%
%  where [a,b]x[c,d] is the domain of f. 
% 
% See also CUMSUM, SUM, SUM2.

if ( isempty( f ) ) % check for empty chebfun2.
    f = [];
    return;
end

f.cols = cumsum( f.cols );   % cumsum along the columns.
f.rows = cumsum( f.rows );   % cumsum along the rows.

end