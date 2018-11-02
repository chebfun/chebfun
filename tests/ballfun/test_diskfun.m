function pass = test_diskfun( )

eps = 1e-10;

% Example 1
f = ballfun(@(x,y,z)cos(x.*y),'cart');
g = diskfun(f,'x','y');
h = diskfun(@(x,y)cos(x.*y));
pass(1) = norm(g-h) < eps;

if (nargout > 0)
    pass = all(pass(:));
end
end
