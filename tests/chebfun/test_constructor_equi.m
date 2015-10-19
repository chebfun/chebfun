function pass = test_constructor_equi(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xx = 2 * rand(100, 1) - 1;

% ATs original quasimatrix case, the following should return a quasimatrix 
% of size [-1,1] x 10 with constant-valued columns:
f = @(x) cos(x); 
x = linspace(-1,1,10); 
v = f(x); 
g = chebfun( v, 'equi' );
pass(1) = isequal(size(g), [Inf 10]) && ...
    norm(g(xx,:) - repmat(v, 100, 1)) < 10*vscale(g)*eps;

% Equi should work with short columns:
g = chebfun( v(1:3)', 'equi' );
pass(2) = isequal(size(g), [Inf,1]);

% Equi should work with shorter columns:
g = chebfun( v(1:2)', 'equi' );
pass(3) = isequal(size(g), [Inf,1]);

% Equi should work with matrices:
g = chebfun( repmat(v',1,10), 'equi' );
pass(4) = isequal(size(g), [Inf,10]);

% Equi should work and matrices that are wider than tall:
g = chebfun( repmat(v',1,11), 'equi' );
pass(5) = isequal(size(g), [Inf,11]);

% Equi should make lines:
g = chebfun( [-1 0 1]', 'equi' );
pass(6) = norm(g - chebfun('x')) < 10*vscale(g)*eps;

% Equi should make lines:
g = chebfun([-3 -1 1 3]','equi');
pass(7) = norm(g - chebfun('3*x')) < 10*vscale(g)*eps;

% Equi should make lines:
g = chebfun( [-1e5 1e5]', 'equi' );
pass(8) = norm(g - chebfun('1e5*x')) < 10*vscale(g)*eps;

% Equi should make parabolas:
g = chebfun( [0 1 0]', 'equi' );
pass(9) = norm(g - chebfun('1 - x.^2')) < 10*vscale(g)*eps;

% Equi should also work for non-quasimatrix constant valued functions:
u = xx(end);
g = chebfun( u, 'equi');
pass(10) = norm(g - u) < 10*vscale(g)*eps;

end
