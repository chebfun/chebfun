function pass = test_vscale( ) 

tol = 1e-2;

% Example 1 : x
f = ballfun(@(x,y,z)x);
v = vscale(f);
pass(1) = norm(v - 1) < tol;

% Example 2 : -x^2
f = ballfun(@(x,y,z)-x.^2);
v = vscale(f);
pass(2) = norm(v - 1) < tol;

% Example 3 : 5*y
f = ballfun(@(x,y,z)5*y);
v = vscale(f);
pass(3) = norm(v - 5) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
