function G = potential( f )
%POTENTIAL  2D vector potential of a CHEBFUN2.
%   G = POTENTIAL(F) where F is a CHEBFUN2 returns a vector-valued CHEBFUN2V
%   with two components such that F = curl(G).
% 
%   Note this is NOT the 3D vector potential because CHEBFUN2 represents
%   functions with two variables.
%
%   TODO: This function is slow and requires improvements. It works for small
%   degree bivariate polynomials.
% 
% See also CHEBFUN2V/CURL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) )
    G = chebfun2v( ); 
    return
end

dom = f.domain; 
x = chebfun2(@(x,y) x, dom); 
y = chebfun2(@(x,y) y, dom); 

% One can show that:
%      f(x,y) = dQ/dx - dP/dy, 
% where Q = x.*S(x,y), P = -y.*S(x,y), and 
%      S(x,y) = integral( s.*f(s.*x,s.*y), [0 1] ) 

S = @(x, y) sum( chebfun(@(s) feval(f, s.*x, s.*y).*s, [0 1] ));
S = chebfun2(S, dom, 'vectorize');

Q = x.*S; 
P = -y.*S; 
G = [P ; Q];

end
