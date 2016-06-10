function pass = test_ndf() 
% Test CHEBFUN3T/NDF

j = 1;

f = chebfun3t();
pass(j) = ndf(f) == 0;
j = j+1;

ff = @(x,y,z) 10;
f = chebfun3t(ff);
pass(j) = ndf(f) == 1;
j = j+1;

f = chebfun3t(@(x,y,z) x, [-1 2 -pi/2 pi -3 1]); 
pass(j) = ndf(f) == 2;
j = j+1;

f = chebfun3t(@(x,y,z) y, [-1 2 -pi/2 pi -3 1]); 
pass(j) = ndf(f) == 2;
j = j+1;

f = chebfun3t(@(x,y,z) z, [-1 2 -pi/2 pi -3 1]); 
pass(j) = ndf(f) == 2;
j = j+1;

ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3t(ff);
ndf_f = prod(size(f.coeffs));
pass(j) = ndf(f) == ndf_f;
j = j+1;

end