function pass = test_ndf() 
% Test CHEBFUN3T/NDF

f = chebfun3t();
pass(1) = ndf(f) == 0;

ff = @(x,y,z) 10;
f = chebfun3t(ff);
pass(2) = ndf(f) == 1;

f = chebfun3t(@(x,y,z) x, [-1 2 -pi/2 pi -3 1]); 
pass(3) = ndf(f) == 2;

f = chebfun3t(@(x,y,z) y, [-1 2 -pi/2 pi -3 1]); 
pass(4) = ndf(f) == 2;

f = chebfun3t(@(x,y,z) z, [-1 2 -pi/2 pi -3 1]); 
pass(5) = ndf(f) == 2;

ff = @(x,y,z) sin(pi*(x+y+z));
f = chebfun3t(ff);
ndf_f = prod(size(f.coeffs));
pass(6) = ndf(f) == ndf_f;

end