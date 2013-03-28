% Test file for chebtech2/chebpts.m

function pass = test_chebpts(varargin)

% Set a tolerance (pref.eps doesn't matter)
tol = 10*eps;

% Test that n = 0 returns empty results:
[x, w, v] = chebtech2.chebpts(0);
pass(1) = ( isempty(x) && isempty(w) && isempty(v) );

% Test n = 1:
[x, w, v] = chebtech2.chebpts(1);
pass(2) = ( x == 0 && w == 2 && v == 1);

% Test that n = 2 returns [-1 ; 1]:
x = chebtech2.chebpts(2);
pass(3) = ( all(size(x) == [2, 1]) && all(x == [-1 ; 1]) );
[x, w, v] = chebtech2.chebpts(2);
pass(4) = ( all(size(w) == [1, 2]) && all(w == [1, 1]) );
pass(5) = ( all(size(v) == [2, 1]) && all(v == .5*[-1 ; 1]) );

% Test that n = 3 returns [-1 ; 0 ; 1]:
x = chebtech2.chebpts(3);
pass(6) = ( all(size(x) == [3, 1]) && all(x == [-1 ; 0 ; 1]) );
[x, w, v] = chebtech2.chebpts(3);
pass(7) = ( all(size(w) == [1, 3]) && norm( w - ([0, 1, 0] + 1/3), inf) < tol );
pass(8) = ( all(size(v) == [3, 1]) && norm( v - ([.5 ; -1 ; .5]) , inf) < tol );

% Test that n = 129 returns vectors of the correct size:
n = 129;
[x, w, v] = chebtech2.chebpts(n);
pass(9) = ( all([size(x) == [n, 1], size(w) == [1, n], size(v) == [n, 1]]) );
% and that the nodes are symmetric:
pass(10) = ( norm(x(1:(n-1)/2) + x(n:-1:(n+3)/2), inf) == 0 );
pass(11) =  ( x((n+1)/2) == 0 );

end
