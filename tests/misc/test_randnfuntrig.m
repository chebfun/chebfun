function pass = test_randnfuntrig( pref )

if ( nargin == 0 ) 
    pref = chebfunpref();
end

rng(0)
f = randnfuntrig(.01);
pass(1) = (abs(std(f)-1) < .1);
pass(2) = (abs(mean(f)) < .1);

rng(0)
g = randnfuntrig(.01,1,[5,7]);
pass(3) = ((sum(f)-sum(g)) < 1e-3);

A = randnfuntrig(1,10,[0 10]);
X = cov(A);
pass(4) = (norm(X-X') == 0);

f = randnfuntrig(1) + 1i*randnfuntrig(1);
pass(5) = (abs(f(-.999)-f(.999)) < .1);

rng(0), f1 = randnfuntrig(1/16,'norm')/4;
rng(0), f2 = randnfuntrig(1/16);
pass(6) = norm(f1-f2) == 0;

end
