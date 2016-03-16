% Test file for chebtech1/chebpts.m

function pass = test_chebpts(varargin)

% Set a tolerance (pref.chebfuneps doesn't matter)
tol = 10*eps;

% Test that n = 0 returns empty results:
[x, w, v] = chebtech1.chebpts(0);
pass(1) = ( isempty(x) && isempty(w) && isempty(v) );

% Test n = 1:
[x, w, v] = chebtech1.chebpts(1);
pass(2) = ( x == 0 && w == 2 && v == 1);

% Test that when n = 2, x = v = [-1/sqrt(2) ; 1/sqrt(2)] and w = :
x = chebtech1.chebpts(2);
pass(3) = ( all(size(x) == [2, 1]) && norm( x - [-1/sqrt(2) ; 1/sqrt(2)], inf) < tol );
[ignored, w] = chebtech1.chebpts(2);
pass(4) = ( all(size(w) == [1, 2]) && all( w == [1, 1]) );
[ignored1, ignored2, v] = chebtech1.chebpts(2);
pass(5) = ( all(size(v) == [2, 1]) && norm( v - [-1/sqrt(2) ; 1/sqrt(2)], inf) < tol );

% Test that n = 3 returns [-1 ; 0 ; 1]:
x = chebtech1.chebpts(3);
pass(6) = ( all(size(x) == [3, 1]) && norm(x - [-sqrt(3)/2 ; 0 ; sqrt(3)/2], inf) < tol );
[ignored, w] = chebtech1.chebpts(3);
pass(7) = ( all(size(w) == [1, 3]) && norm( w - ([4/9, 10/9, 4/9]), inf) < tol );
[ignored1, ignored2, v] = chebtech1.chebpts(3);
pass(8) = ( all(size(v) == [3, 1]) && norm( v - ([.5 ; -1 ; .5]) , inf) < tol );

% Test that n = 129 returns vectors of the correct size:
n = 129;
[x, w, v] = chebtech1.chebpts(n);
pass(9) = ( all([size(x) == [n, 1], size(w) == [1, n], size(v) == [n, 1]]) );
% and that the nodes are symmetric:
pass(10) = ( norm(x(1:(n-1)/2) + x(n:-1:(n+3)/2), inf) < tol );
pass(11) =  ( x((n+1)/2) < tol );

end
