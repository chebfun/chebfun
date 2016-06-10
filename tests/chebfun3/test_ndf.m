function pass = test_ndf() 
% Test CHEBFUN3/NDF

j = 1;

f = chebfun3();
pass(j) = ndf(f) == 0;
j = j+1;

ff = @(x,y,z) 10;
f = chebfun3(ff);
pass(j) = ndf(f) == 4;
j = j+1;

f = chebfun3(@(x,y,z) x, [-1 2 -pi/2 pi -3 1]); 
pass(j) = ndf(f) == 5;
j = j+1;

f = chebfun3(@(x,y,z) y); 
pass(j) = ndf(f) == 5;
j = j+1;

f = chebfun3(@(x,y,z) z, [-1 2 -pi/2 pi -3 1]); 
pass(j) = ndf(f) == 5;
j = j+1;

ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3(ff);
[r1, r2, r3] = rank(f);
[m, n, p] = length(f);
ndf_f = r1*m + r2*n + r3*p + r1*r2*r3;
pass(j) = ndf(f) == ndf_f;
j = j+1;

end