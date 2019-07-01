function pass = test_end(pref)

if ( nargin == 0) 
    pref = chebfunpref; 
end

% See #2282.
u = chebfun2(@(x, y) x+y, [0 2 0 3], pref);

pass(1) = abs(u(end,end) - u(2,3)) < 1e-14;

x = chebfun('x', [0 2]);
f = x' + 3;
pass(2) = norm(u(:,end) - f) < 1e-14;

y = chebfun('x', [0 3]);
g = 2 + y;
pass(3) = norm(u(end,:) - g) < 1e-14;

end