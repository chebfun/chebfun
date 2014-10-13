% Test file for @chebfun/polyval.m.

function pass = test_polyval(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(666)

%% Test scalar problem:

c = rand(3,1);
x = chebfun('x', [0 pi], pref);
f = polyval(c, x);
err = norm(c(1)*x.^2 + c(2)*x + c(3) - f, inf);
tol = vscale(f).*epslevel(f);
pass(1) = err < 10*tol;

%% Test quasimatrix:

x1 = x;
x2 = x.^2;
xx = [x1, x2];
f2 = polyval(c, xx);
g1 = c(1)*x1.^2 + c(2)*x1 + c(3);
g2 = c(1)*x2.^2 + c(2)*x2 + c(3);
err1 = norm(g1 - f2(:,1), inf);
err2 = norm(g2 - f2(:,2), inf);
pass(2) = all([err1, err2] < 10*epslevel(f2).*vscale(f2));

%% Test matrix:

seedRNG(667)
c2 = rand(3,1);
cc = [c, c2];
f3 = polyval(cc, x);
g4 = c2(1)*x.^2 + c2(2)*x + c2(3);
err3 = norm(g1 - f3(:,1), inf);
err4 = norm(g4 - f3(:,2), inf);
pass(3) = all([err3, err4] < 10*epslevel(f3).*vscale(f3));

%% Test matrix and quasimatrix:

try 
    polyval(cc, xx);
    pass(4) = 0;
catch
    pass(4) = true;
end

end
