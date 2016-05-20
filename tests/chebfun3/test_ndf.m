function pass = test_ndf() 
% Test CHEBFUN3/NDF

f = chebfun3();
pass(1) = ndf(f) == 0;

ff = @(x,y,z) 10;
f = chebfun3(ff);
pass(2) = ndf(f) == 4;

f = chebfun3(@(x,y,z) x, [-1 2 -pi/2 pi -3 1]); 
pass(3) = ndf(f) == 5;

f = chebfun3(@(x,y,z) y); 
pass(4) = ndf(f) == 5;

f = chebfun3(@(x,y,z) z, [-1 2 -pi/2 pi -3 1]); 
pass(5) = ndf(f) == 5;

ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff);
[rX, rY, rZ] = rank(f);
[m, n, p] = length(f);
ndf_f = rX*m + rY*n + rZ*p + rX*rY*rZ;
pass(6) = ndf(f) == ndf_f;

end