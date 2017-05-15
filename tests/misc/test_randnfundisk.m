function pass = test_randnfundisk( pref )

if ( nargin == 0 ) 
    pref = chebfunpref();
end

rng(0)
f = randnfundisk(.2);
pass(1) = (abs(mean2(f.^2)-1) < .1);
pass(2) = (abs(mean2(f)) < .1);

f = randnfundisk(40);
d = max2(f) - min2(f);
pass(3) = (d < 0.5) & (d > 0.001);



