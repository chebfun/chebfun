function pass = test_randnfunsphere( pref )

if ( nargin == 0 ) 
    pref = chebfunpref();
end

rng(0)
f = randnfunsphere(.2);
pass(1) = (abs(mean2(f.^2)-1) < .1);
pass(2) = (abs(mean2(f)) < .1);

f = randnfunsphere(1e6);
pass(3) = norm(diff(f),'fro') < 1e-4;

f = randnfunsphere(3.1);
pass(4) = ( rank(f) == 4 );

f = randnfunsphere(3.1,'mono');
pass(5) = ( rank(f) == 3 );
